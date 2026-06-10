"""
sota_models.py
==============
Faithful re-implementations of recent sequence-based antibody-antigen affinity SOTA, to be trained
on the SAME SAaIntDB data, SAME 10-fold splits and SAME frozen ESM-2 650M features as our model
(none of these methods released public code, so we reimplement from the papers).

* DuaDeepSeqAffinity  -- Dual-stream CNN + Transformer per protein (arXiv:2512.22007).
    paper config: per-residue ESM-2 650M (1280); per protein a Transformer branch (GAP->1280) and a
    CNN branch [Conv1d 256/k3 -> Conv1d 128/k5 -> GAP=128]; intra-fuse concat (1408); inter-fuse
    concat antibody+antigen (2816) -> FC head -> 1; MSE.
* DGAffinityConvNeXt  -- PLM features + 1D-ConvNeXt stack per protein -> GAP -> concat -> FC -> 1
    (DG-Affinity family; arch approximated from the paper, ConvNeXt-style depthwise blocks).

`reduced=True` projects 1280->d_model (default 256) and uses a 1-layer Transformer so the models are
trainable on CPU; `reduced=False` uses the paper dimensions (recommended on GPU). The choice is logged.
"""
import torch, torch.nn as nn, torch.nn.functional as F


def masked_mean(x, mask):
    # x:(B,L,C)  mask:(B,L) 1=valid
    m = mask.unsqueeze(-1).float()
    return (x * m).sum(1) / m.sum(1).clamp(min=1.0)


class ProteinDualStream(nn.Module):
    """Transformer branch + CNN branch over per-residue embeddings -> [GAP_tr ; GAP_cnn]."""
    def __init__(self, in_dim=1280, reduced=True, d_model=256, n_heads=4, n_layers=1):
        super().__init__()
        self.reduced = reduced
        dm = d_model if reduced else in_dim
        self.proj = nn.Linear(in_dim, dm) if reduced else nn.Identity()
        enc = nn.TransformerEncoderLayer(d_model=dm, nhead=n_heads,
                                         dim_feedforward=4 * dm, batch_first=True, dropout=0.1)
        self.tr = nn.TransformerEncoder(enc, num_layers=n_layers)
        self.c1 = nn.Conv1d(dm, 256, kernel_size=3, padding=1)
        self.c2 = nn.Conv1d(256, 128, kernel_size=5, padding=2)
        self.out_dim = dm + 128

    def forward(self, x, mask):                       # x:(B,L,1280) mask:(B,L)
        x = self.proj(x)
        pad = ~mask.bool()
        t = self.tr(x, src_key_padding_mask=pad)      # (B,L,dm)
        t = masked_mean(t, mask)                      # (B,dm)
        c = F.relu(self.c1(x.transpose(1, 2)))
        c = F.relu(self.c2(c)).transpose(1, 2)        # (B,L,128)
        c = masked_mean(c, mask)                      # (B,128)
        return torch.cat([t, c], dim=-1)


class DuaDeepSeqAffinity(nn.Module):
    def __init__(self, reduced=True, d_model=256):
        super().__init__()
        self.ab = ProteinDualStream(reduced=reduced, d_model=d_model)
        self.ag = ProteinDualStream(reduced=reduced, d_model=d_model)
        d = self.ab.out_dim + self.ag.out_dim
        self.head = nn.Sequential(
            nn.Linear(d, 512), nn.ReLU(), nn.Dropout(0.2),
            nn.Linear(512, 128), nn.ReLU(), nn.Dropout(0.2),
            nn.Linear(128, 1))

    def forward(self, ab_x, ab_m, ag_x, ag_m):
        f = torch.cat([self.ab(ab_x, ab_m), self.ag(ag_x, ag_m)], dim=-1)
        return self.head(f).squeeze(-1)


class ConvNeXt1DBlock(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.dw = nn.Conv1d(dim, dim, kernel_size=7, padding=3, groups=dim)
        self.norm = nn.LayerNorm(dim)
        self.pw1 = nn.Linear(dim, 4 * dim); self.pw2 = nn.Linear(4 * dim, dim)

    def forward(self, x):                              # x:(B,L,C)
        h = self.dw(x.transpose(1, 2)).transpose(1, 2)
        h = self.norm(h); h = self.pw2(F.gelu(self.pw1(h)))
        return x + h


class DGAffinityConvNeXt(nn.Module):
    """PLM features -> per-protein 1D-ConvNeXt stack -> masked GAP -> concat -> FC -> 1."""
    def __init__(self, in_dim=1280, reduced=True, d_model=256, depth=3):
        super().__init__()
        dm = d_model if reduced else 512
        self.proj = nn.Linear(in_dim, dm)
        self.ab = nn.ModuleList([ConvNeXt1DBlock(dm) for _ in range(depth)])
        self.ag = nn.ModuleList([ConvNeXt1DBlock(dm) for _ in range(depth)])
        self.head = nn.Sequential(nn.Linear(2 * dm, 256), nn.ReLU(), nn.Dropout(0.2),
                                  nn.Linear(256, 1))

    def _enc(self, x, mask, blocks):
        x = self.proj(x)
        for b in blocks:
            x = b(x)
        return masked_mean(x, mask)

    def forward(self, ab_x, ab_m, ag_x, ag_m):
        f = torch.cat([self._enc(ab_x, ab_m, self.ab), self._enc(ag_x, ag_m, self.ag)], dim=-1)
        return self.head(f).squeeze(-1)


MODELS = {'DuaDeep-SeqAffinity': DuaDeepSeqAffinity, 'DG-Affinity': DGAffinityConvNeXt}
