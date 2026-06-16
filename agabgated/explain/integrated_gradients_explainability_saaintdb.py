"""
Integrated Gradients attribution for MutualTriStreamStrong — SAaIntDB dataset.
Per-residue ESM-2 embeddings → IG attribution → bar plots + rectangular heatmaps.

Sequences are loaded directly from the SAaIntDB CSV (no FASTA needed).
Nanobodies (empty L_seq) reuse the heavy-chain embedding for the light slot,
matching the training-time convention used in mutual_strong_saaintdb.py.

Figures per complex:
  1. Bar plots (per chain)       → bar_{complex}_{chain}.png/pdf
  2. Per-chain rectangular maps  → heatmap_{complex}_{chain}.png/pdf
  3. Combined 3-panel figure     → combined_{complex}.png/pdf
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
from captum.attr import IntegratedGradients
from transformers import AutoModel, AutoTokenizer

from agabgated.models.mutual_strong import MutualTriStreamStrong
from agabgated.utils.main_symmetric_mean import cosine_to_affinity

# =============================================================================
# Config
# =============================================================================

DEVICE     = 'cuda' if torch.cuda.is_available() else 'cpu'
MODEL_PATH = "model_weights/saaintdb_allcdr_random_bestfold.pt"
CSV_PATH   = "data/saaintdb_with_antigen_names.csv"
ESM2_MODEL = "facebook/esm2_t33_650M_UR50D"
IG_STEPS   = 50
OUTPUT_DIR = "results_mutual_strong_saaintdb/ig_attribution"
N_COMPLEXES = 5           # how many complexes to sample from SAaIntDB
os.makedirs(OUTPUT_DIR, exist_ok=True)

_INVALID_SEQS = {'', 'nan', 'N/A', 'NA', 'N.A.', 'None'}

def _is_valid(val) -> bool:
    if pd.isna(val):
        return False
    return str(val).strip() not in _INVALID_SEQS


# =============================================================================
# Load representative complexes from SAaIntDB CSV
# =============================================================================

def load_complexes(csv_path: str, n: int = N_COMPLEXES) -> dict:
    """
    Pick N unique PDB_IDs from the SAaIntDB CSV and return their sequences.
    For nanobodies (empty L_seq), the heavy sequence is reused in the light slot,
    matching what mutual_strong_saaintdb.py does during training.
    """
    df = pd.read_csv(csv_path)
    # Filter: must have H_seq and Ag_seq
    df = df[df['H_seq'].apply(_is_valid) & df['Ag_seq'].apply(_is_valid)]
    # One row per PDB_ID (first occurrence)
    df_unique = df.drop_duplicates('PDB_ID').head(n)

    complexes = {}
    for _, row in df_unique.iterrows():
        pdb = str(row['PDB_ID'])
        heavy_seq = str(row['H_seq']).strip()
        ag_seq    = str(row['Ag_seq']).strip()

        if _is_valid(row.get('L_seq', '')):
            light_seq = str(row['L_seq']).strip()
            is_nanobody = False
        else:
            light_seq   = heavy_seq   # nanobody: reuse heavy
            is_nanobody = True

        complexes[pdb] = {
            'heavy':   heavy_seq,
            'light':   light_seq,
            'antigen': ag_seq,
            'nanobody': is_nanobody,
        }
        print(f"  {pdb}: H={len(heavy_seq)} L={len(light_seq)} "
              f"Ag={len(ag_seq)} {'[nanobody]' if is_nanobody else ''}")
    return complexes


# =============================================================================
# ESM-2 per-residue embedder
# =============================================================================

class ESM2PerResidueEmbedder:
    """Load ESM-2 once; return per-residue embeddings (CLS/EOS stripped)."""

    def __init__(self, model_name=ESM2_MODEL, device=DEVICE):
        print(f"  Loading ESM-2 ({model_name}) ...")
        self.device    = device
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.esm       = AutoModel.from_pretrained(model_name).to(device)
        self.esm.eval()
        print("  ESM-2 ready.")

    @torch.no_grad()
    def embed(self, sequence: str) -> np.ndarray:
        """Return (L, 1280) float32 array — CLS and EOS tokens removed."""
        tokens = self.tokenizer(
            sequence, return_tensors='pt',
            padding=False, truncation=True, max_length=1024
        ).to(self.device)
        out = self.esm(**tokens)
        emb = out.last_hidden_state.squeeze(0)[1:-1].cpu().float().numpy()
        return emb

    def cleanup(self):
        del self.esm, self.tokenizer
        torch.cuda.empty_cache()


# =============================================================================
# Per-residue model wrapper for IG
# =============================================================================

class MutualTriStreamPerResidue(torch.nn.Module):
    """
    Wraps MutualTriStreamStrong for per-residue attribution.
    Inputs : light_emb (B, L_l, D), heavy_emb (B, L_h, D), antigen_emb (B, L_a, D)
    Output : predicted pKd as scalar
    """
    def __init__(self, model: MutualTriStreamStrong, pkd_bounds):
        super().__init__()
        self.model  = model
        self.pkd_lo = pkd_bounds[0]
        self.pkd_hi = pkd_bounds[1]

    def forward(self, light_emb, heavy_emb, antigen_emb):
        le  = light_emb.mean(dim=1)
        he  = heavy_emb.mean(dim=1)
        ae  = antigen_emb.mean(dim=1)
        out = self.model(le, he, ae)
        cos = out['cosine_similarity']
        pkd = (cos + 1.0) / 2.0 * (self.pkd_hi - self.pkd_lo) + self.pkd_lo
        return pkd.squeeze()


# =============================================================================
# Model loading
# =============================================================================

def load_model(model_path: str, device=DEVICE):
    ckpt   = torch.load(model_path, map_location=device)
    cfg    = ckpt.get('config', {})
    model  = MutualTriStreamStrong(
        esm_dim        = cfg.get('esm_dim', 1280),
        projected_size = cfg.get('projected_size', 256),
        num_heads      = cfg.get('num_heads', 8),
        dropout        = 0.0,
        n_layers       = cfg.get('n_layers', 1),
        device         = device
    ).to(device)
    model.load_state_dict(ckpt['model_state_dict'])
    model.eval()
    # SAaIntDB checkpoints store bounds under 'pkd_bounds'
    pkd_bounds = ckpt.get('pkd_bounds', ckpt.get('train_bounds', (5.0, 12.0)))
    return model, pkd_bounds


# =============================================================================
# Predict affinity from mean-pooled embeddings
# =============================================================================

@torch.no_grad()
def predict_affinity(model, light_np, heavy_np, antigen_np, pkd_bounds, device=DEVICE):
    le = torch.tensor(light_np.mean(0, keepdims=True),   dtype=torch.float32).to(device)
    he = torch.tensor(heavy_np.mean(0, keepdims=True),   dtype=torch.float32).to(device)
    ae = torch.tensor(antigen_np.mean(0, keepdims=True), dtype=torch.float32).to(device)
    out = model(le, he, ae)
    cos = out['cosine_similarity'].cpu().numpy()[0]
    return cosine_to_affinity(cos, *pkd_bounds)


# =============================================================================
# Integrated Gradients
# =============================================================================

def compute_ig_attributions(wrapper, light_np, heavy_np, antigen_np,
                             device=DEVICE, n_steps=IG_STEPS):
    """
    Run IG simultaneously over all three inputs (zero baseline).
    Returns (light_attr, heavy_attr, antigen_attr) each shape (L,) [summed over D].
    """
    le = torch.tensor(light_np[np.newaxis],   dtype=torch.float32,
                      requires_grad=True).to(device)
    he = torch.tensor(heavy_np[np.newaxis],   dtype=torch.float32,
                      requires_grad=True).to(device)
    ae = torch.tensor(antigen_np[np.newaxis], dtype=torch.float32,
                      requires_grad=True).to(device)

    bl_le = torch.zeros_like(le)
    bl_he = torch.zeros_like(he)
    bl_ae = torch.zeros_like(ae)

    ig    = IntegratedGradients(wrapper)
    attrs = ig.attribute(
        inputs    = (le, he, ae),
        baselines = (bl_le, bl_he, bl_ae),
        n_steps   = n_steps,
        return_convergence_delta=False
    )
    light_attr   = attrs[0].squeeze(0).sum(dim=-1).detach().cpu().numpy()
    heavy_attr   = attrs[1].squeeze(0).sum(dim=-1).detach().cpu().numpy()
    antigen_attr = attrs[2].squeeze(0).sum(dim=-1).detach().cpu().numpy()
    return light_attr, heavy_attr, antigen_attr


# =============================================================================
# Plotting helpers
# =============================================================================

def _norm(attr):
    lo, hi = attr.min(), attr.max()
    return (attr - lo) / (hi - lo + 1e-8)


def plot_bar(attr, seq, chain_name, complex_name):
    attr_n = _norm(attr)
    pos    = np.arange(len(attr_n)) + 1

    fig, ax = plt.subplots(figsize=(max(10, len(pos) * 0.12), 4))
    colours = ['#d73027' if v >= 0.6 else '#4575b4' if v <= 0.3 else '#74add1'
               for v in attr_n]
    ax.bar(pos, attr_n, width=0.85, color=colours, alpha=0.85)
    ax.set_xlabel('Residue position', fontsize=10)
    ax.set_ylabel('Normalised IG attribution', fontsize=10)
    ax.set_title(f'{chain_name} chain – {complex_name} (SAaIntDB / MutualTriStreamStrong)',
                 fontsize=11)
    ax.set_xlim(0.5, len(pos) + 0.5)
    ax.set_ylim(0, 1.05)
    ax.axhline(0.5, color='gray', lw=0.7, ls='--', alpha=0.6)
    ax.grid(axis='y', alpha=0.25)
    if len(pos) > 40:
        step = max(1, len(pos) // 20)
        ax.set_xticks(pos[::step])
        ax.set_xticklabels(pos[::step], fontsize=8)
    else:
        ax.set_xticks(pos)
        ax.set_xticklabels(list(seq), fontsize=7)

    import matplotlib.patches as mpatches
    legend_patches = [
        mpatches.Patch(color='#d73027', label='High (>0.6)'),
        mpatches.Patch(color='#74add1', label='Mid'),
        mpatches.Patch(color='#4575b4', label='Low (<0.3)'),
    ]
    ax.legend(handles=legend_patches, fontsize=8, loc='upper right')
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f'bar_{complex_name}_{chain_name}')
    plt.savefig(fname + '.png', dpi=300)
    plt.savefig(fname + '.pdf')
    plt.close(fig)


def plot_rectangular_heatmap(attr, seq, chain_name, complex_name, cols=20):
    cmap_dict = {'Light': 'Blues', 'Heavy': 'Greens', 'Antigen': 'Purples'}
    cmap      = cmap_dict.get(chain_name, 'RdYlGn')

    attr_n = _norm(attr)
    L      = len(attr_n)
    rows   = int(np.ceil(L / cols))
    pad    = rows * cols - L
    p_attr = np.pad(attr_n, (0, pad), constant_values=np.nan)
    p_seq  = list(seq) + [''] * pad
    hm     = p_attr.reshape(rows, cols)
    sg     = np.array(p_seq).reshape(rows, cols)

    fig, ax = plt.subplots(figsize=(cols * 0.45, rows * 0.45 + 0.8))
    im = ax.imshow(hm, cmap=cmap, vmin=0, vmax=1, aspect='auto')

    ax.set_xticks(np.arange(-0.5, cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, rows, 1), minor=True)
    ax.grid(which='minor', color='white', linestyle='-', linewidth=0.6)
    ax.tick_params(which='minor', bottom=False, left=False)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    for i in range(rows):
        for j in range(cols):
            aa = sg[i, j]
            if aa != '':
                txt_col = 'white' if hm[i, j] > 0.6 else 'black'
                ax.text(j, i, aa, ha='center', va='center',
                        color=txt_col, fontsize=6.5, fontweight='bold')

    ax.set_title(f'{chain_name} chain – {complex_name} (SAaIntDB)', fontsize=10, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cbar.set_label('Normalised IG importance', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f'heatmap_{complex_name}_{chain_name}')
    plt.savefig(fname + '.png', dpi=300)
    plt.savefig(fname + '.pdf')
    plt.close(fig)


def plot_combined_heatmaps(attrs, seqs, chain_names, complex_name,
                            cols=20, pred_pkd=None, nanobody=False):
    cmaps = ['Blues', 'Greens', 'Purples']
    n     = len(attrs)
    fig, axes = plt.subplots(n, 1, figsize=(max(16, cols * 0.5), n * 3.5))
    if n == 1:
        axes = [axes]
    fig.patch.set_facecolor('#f8f8f8')
    plt.subplots_adjust(hspace=0.5)

    last_im = None
    for idx, (attr, seq, chain_name, ax) in enumerate(
            zip(attrs, seqs, chain_names, axes)):
        attr_n = _norm(attr)
        L      = len(attr_n)
        rows   = int(np.ceil(L / cols))
        pad    = rows * cols - L
        p_attr = np.pad(attr_n, (0, pad), constant_values=np.nan)
        p_seq  = list(seq) + [''] * pad
        hm     = p_attr.reshape(rows, cols)
        sg     = np.array(p_seq).reshape(rows, cols)

        last_im = ax.imshow(hm, cmap=cmaps[idx], vmin=0, vmax=1, aspect='auto')
        ax.set_xticks(np.arange(-0.5, cols, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, rows, 1), minor=True)
        ax.grid(which='minor', color='white', linestyle='-', linewidth=0.5)
        ax.tick_params(which='minor', bottom=False, left=False)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        for i in range(rows):
            for j in range(cols):
                aa = sg[i, j]
                if aa != '':
                    txt_col = 'white' if hm[i, j] > 0.6 else 'black'
                    ax.text(j, i, aa, ha='center', va='center',
                            color=txt_col, fontsize=6.5, fontweight='bold')
        ax.set_title(f'Residue Importance: {chain_name}', fontsize=12, fontweight='bold')
        ax.set_facecolor('#f0f0f0')

    nanobody_note = ' [Nanobody — heavy reused for light]' if nanobody else ''
    title_extra   = f'  |  Predicted pKd = {pred_pkd:.2f}' if pred_pkd is not None else ''
    fig.suptitle(
        f'Residue Importance (IG) – {complex_name} (SAaIntDB){nanobody_note}{title_extra}',
        fontsize=13, fontweight='bold'
    )
    cbar = fig.colorbar(last_im, ax=axes, fraction=0.018, pad=0.02)
    cbar.set_label('Normalised IG attribution', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fname = os.path.join(OUTPUT_DIR, f'combined_{complex_name}')
    plt.savefig(fname + '.png', dpi=300, bbox_inches='tight')
    plt.savefig(fname + '.pdf',            bbox_inches='tight')
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================

def main():
    print(f"Device: {DEVICE}")

    # ── load complexes ────────────────────────────────────────────────────────
    print(f"\nLoading {N_COMPLEXES} SAaIntDB complexes from {CSV_PATH} ...")
    complexes = load_complexes(CSV_PATH, n=N_COMPLEXES)
    if not complexes:
        print("No valid complexes found. Exiting.")
        return

    # ── load model ────────────────────────────────────────────────────────────
    print(f"\nLoading MutualTriStreamStrong from {MODEL_PATH} ...")
    model, pkd_bounds = load_model(MODEL_PATH)
    print(f"  pkd_bounds: {pkd_bounds}")

    wrapper = MutualTriStreamPerResidue(model, pkd_bounds).to(DEVICE)
    wrapper.eval()

    # ── load ESM-2 ────────────────────────────────────────────────────────────
    embedder = ESM2PerResidueEmbedder(ESM2_MODEL, DEVICE)

    # ── process each complex ──────────────────────────────────────────────────
    for name, info in complexes.items():
        print(f"\n{'='*60}")
        print(f"Complex: {name}  {'[nanobody]' if info['nanobody'] else ''}")
        print(f"{'='*60}")

        print("  Embedding sequences ...")
        light_emb   = embedder.embed(info['light'])    # (L_l, 1280)
        heavy_emb   = embedder.embed(info['heavy'])    # (L_h, 1280)
        antigen_emb = embedder.embed(info['antigen'])  # (L_a, 1280)
        print(f"  Light: {light_emb.shape[0]} res | "
              f"Heavy: {heavy_emb.shape[0]} res | "
              f"Antigen: {antigen_emb.shape[0]} res")

        pred = predict_affinity(model, light_emb, heavy_emb, antigen_emb, pkd_bounds)
        print(f"  Predicted pKd: {pred:.3f}")

        print(f"  Running IG ({IG_STEPS} steps) ...")
        light_attr, heavy_attr, antigen_attr = compute_ig_attributions(
            wrapper, light_emb, heavy_emb, antigen_emb
        )
        print(f"  Attribution ranges: "
              f"Light [{light_attr.min():.3f}, {light_attr.max():.3f}]  "
              f"Heavy [{heavy_attr.min():.3f}, {heavy_attr.max():.3f}]  "
              f"Antigen [{antigen_attr.min():.3f}, {antigen_attr.max():.3f}]")

        print("  Generating figures ...")
        chain_label = 'Light' if not info['nanobody'] else 'Light(=Heavy)'

        plot_bar(light_attr,   info['light'],   chain_label, name)
        plot_bar(heavy_attr,   info['heavy'],   'Heavy',     name)
        plot_bar(antigen_attr, info['antigen'], 'Antigen',   name)

        plot_rectangular_heatmap(light_attr,   info['light'],   chain_label, name)
        plot_rectangular_heatmap(heavy_attr,   info['heavy'],   'Heavy',     name)
        plot_rectangular_heatmap(antigen_attr, info['antigen'], 'Antigen',   name)

        plot_combined_heatmaps(
            [light_attr, heavy_attr, antigen_attr],
            [info['light'], info['heavy'], info['antigen']],
            [chain_label, 'Heavy chain', 'Antigen'],
            name,
            pred_pkd=pred,
            nanobody=info['nanobody']
        )
        print(f"  Saved figures for {name} → {OUTPUT_DIR}/")

    embedder.cleanup()
    print(f"\nDone. All outputs in: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
