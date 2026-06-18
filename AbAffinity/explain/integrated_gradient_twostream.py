"""
Integrated Gradients attribution for TwoStreamSymmetric.
Per-residue ESM-2 embeddings → IG attribution → bar plots + rectangular heatmaps.

The two-stream model treats the antibody as a single concatenated heavy+light
sequence (mean-pooled).  For per-residue attribution we:
  1. Embed heavy and light chains separately via ESM-2.
  2. Concatenate them to form the full antibody residue tensor.
  3. Run IG over (antibody_emb, antigen_emb).
  4. Split the resulting antibody attribution back into heavy and light parts.

Figures per complex:
  1. Bar plots (per chain)       → bar_{complex}_{chain}.png/pdf
  2. Per-chain rectangular maps  → heatmap_{complex}_{chain}.png/pdf
  3. Combined 3-panel figure     → combined_{complex}.png/pdf
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import torch
from captum.attr import IntegratedGradients
from transformers import AutoModel, AutoTokenizer

from AbAffinity.models.two_stream_mutualstrong import TwoStreamSymmetric

# =============================================================================
# Config
# =============================================================================

DEVICE     = 'cuda' if torch.cuda.is_available() else 'cpu'
MODEL_PATH = "models/twostream_balm_seed314.pt"
ESM2_MODEL = "facebook/esm2_t33_650M_UR50D"
IG_STEPS   = 50          # Riemann steps for IG
OUTPUT_DIR = "results_twostream/ig_attribution"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# ESM-2 per-residue embedder  (loaded once)
# =============================================================================

class ESM2PerResidueEmbedder:
    """Load ESM-2 once; return (L, 1280) per-residue embeddings (CLS/EOS stripped)."""

    def __init__(self, model_name=ESM2_MODEL, device=DEVICE):
        print(f"  Loading ESM-2 ({model_name}) …")
        self.device    = device
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.esm       = AutoModel.from_pretrained(model_name).to(device)
        self.esm.eval()
        print("  ESM-2 ready.")

    @torch.no_grad()
    def embed(self, sequence: str) -> np.ndarray:
        """Return (L, 1280) float32 array with CLS/EOS tokens removed."""
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
# Per-residue wrapper for TwoStreamSymmetric
#
# TwoStreamSymmetric.forward(antibody_emb, antigen_emb) expects mean-pooled
# (B, D) tensors.  This wrapper accepts (B, L, D) per-residue tensors, applies
# mean-pooling (so gradients flow through it), and then calls the model.
# =============================================================================

class TwoStreamPerResidue(torch.nn.Module):
    """
    Wraps TwoStreamSymmetric for per-residue IG attribution.
    Inputs : antibody_emb (B, L_ab, D),  antigen_emb (B, L_ag, D)
    Output : predicted pKD as scalar
    """
    def __init__(self, model: TwoStreamSymmetric, pkd_bounds):
        super().__init__()
        self.model  = model
        self.pkd_lo = pkd_bounds[0]
        self.pkd_hi = pkd_bounds[1]

    def forward(self, antibody_emb, antigen_emb):
        ab = antibody_emb.mean(dim=1)   # (B, D)
        ag = antigen_emb.mean(dim=1)    # (B, D)
        out = self.model(ab, ag)
        cos = out['cosine_similarity']
        # Differentiable cosine → pKD conversion
        pkd = (cos + 1.0) / 2.0 * (self.pkd_hi - self.pkd_lo) + self.pkd_lo
        return pkd.squeeze()


# =============================================================================
# Model loading
# =============================================================================

def load_model(model_path: str, device=DEVICE):
    ckpt       = torch.load(model_path, map_location=device)
    cfg        = ckpt.get('config', {})
    model      = TwoStreamSymmetric(
        esm_dim        = 1280,
        projected_size = cfg.get('projected_size', 256),
        num_heads      = cfg.get('num_heads', 8),
        dropout        = 0.0,          # no dropout during attribution
        device         = device
    ).to(device)
    model.load_state_dict(ckpt['model_state_dict'])
    model.eval()
    pkd_bounds = ckpt.get('train_bounds', (5.0, 12.0))
    return model, pkd_bounds


# =============================================================================
# Affinity prediction (reporting only)
# =============================================================================

@torch.no_grad()
def predict_affinity(model, heavy_np, light_np, antigen_np, pkd_bounds, device=DEVICE):
    """Predict pKD from mean-pooled embeddings of concatenated heavy+light and antigen."""
    # Concatenate heavy and light, then mean-pool to match training representation
    ab_np = np.concatenate([heavy_np, light_np], axis=0)
    ab = torch.tensor(ab_np.mean(0, keepdims=True), dtype=torch.float32).to(device)
    ag = torch.tensor(antigen_np.mean(0, keepdims=True), dtype=torch.float32).to(device)
    out = model(ab, ag)
    cos = out['cosine_similarity'].cpu().numpy()[0]
    pkd_lo, pkd_hi = pkd_bounds
    return float((cos + 1.0) / 2.0 * (pkd_hi - pkd_lo) + pkd_lo)


# =============================================================================
# Integrated Gradients (two inputs: antibody, antigen)
# =============================================================================

def compute_ig_attributions(wrapper, heavy_np, light_np, antigen_np,
                             device=DEVICE, n_steps=IG_STEPS):
    """
    Run IG over (antibody_emb, antigen_emb).
    Antibody = concatenation of [heavy | light] residue embeddings.
    Baseline: all-zero embeddings.

    Returns:
        heavy_attr   (L_h,)  – IG attributions for heavy-chain residues
        light_attr   (L_l,)  – IG attributions for light-chain residues
        antigen_attr (L_a,)  – IG attributions for antigen residues
    """
    L_h = heavy_np.shape[0]

    # Concatenate heavy + light to form the full antibody input
    ab_np = np.concatenate([heavy_np, light_np], axis=0)

    ab = torch.tensor(ab_np[np.newaxis],   dtype=torch.float32,
                      requires_grad=True).to(device)   # (1, L_h+L_l, D)
    ag = torch.tensor(antigen_np[np.newaxis], dtype=torch.float32,
                      requires_grad=True).to(device)   # (1, L_a, D)

    bl_ab = torch.zeros_like(ab)
    bl_ag = torch.zeros_like(ag)

    ig    = IntegratedGradients(wrapper)
    attrs = ig.attribute(
        inputs    = (ab, ag),
        baselines = (bl_ab, bl_ag),
        n_steps   = n_steps,
        return_convergence_delta=False
    )

    # attrs[0]: (1, L_h+L_l, D) → sum over D → (L_h+L_l,)
    ab_attr      = attrs[0].squeeze(0).sum(dim=-1).detach().cpu().numpy()
    antigen_attr = attrs[1].squeeze(0).sum(dim=-1).detach().cpu().numpy()

    heavy_attr = ab_attr[:L_h]
    light_attr = ab_attr[L_h:]
    return heavy_attr, light_attr, antigen_attr


# =============================================================================
# Plotting helpers  (identical API to the three-stream script)
# =============================================================================

def _norm(attr):
    """Min-max normalise to [0, 1]."""
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
    ax.set_title(f'{chain_name} chain – {complex_name} (TwoStreamSymmetric)', fontsize=11)
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
        mpatches.Patch(color='#d73027', label='High importance (>0.6)'),
        mpatches.Patch(color='#74add1', label='Mid importance'),
        mpatches.Patch(color='#4575b4', label='Low importance (<0.3)'),
    ]
    ax.legend(handles=legend_patches, fontsize=8, loc='upper right')
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f'bar_{complex_name}_{chain_name}')
    plt.savefig(fname + '.png', dpi=300); plt.savefig(fname + '.pdf')
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

    ax.set_title(f'{chain_name} chain – {complex_name}', fontsize=10, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cbar.set_label('Normalised IG importance', fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, f'heatmap_{complex_name}_{chain_name}')
    plt.savefig(fname + '.png', dpi=300); plt.savefig(fname + '.pdf')
    plt.close(fig)


def plot_combined_heatmaps(attrs, seqs, chain_names, complex_name, cols=20, pred_pkd=None):
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

    title_extra = f'  |  Predicted pKD = {pred_pkd:.2f}' if pred_pkd is not None else ''
    fig.suptitle(f'Residue Importance (IG) – {complex_name}{title_extra}',
                 fontsize=13, fontweight='bold')
    cbar = fig.colorbar(last_im, ax=axes, fraction=0.018, pad=0.02)
    cbar.set_label('Normalised IG attribution', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fname = os.path.join(OUTPUT_DIR, f'combined_{complex_name}')
    plt.savefig(fname + '.png', dpi=300, bbox_inches='tight')
    plt.savefig(fname + '.pdf',            bbox_inches='tight')
    plt.close(fig)


# =============================================================================
# Complexes to analyse  (same set as the three-stream script)
# =============================================================================

COMPLEXES = {
    '1VFB': {
        'light':   "DIVLTQSPASLSASVGETVTITCRASGNIHNYLAWYQQKQGKSPQLLVYYTTTLADGVPSRFSGSGSGTQYSLKINSLQPEDFGSYYCQHFWSTPRTFGGGTKLEIK",
        'heavy':   "QVQLQESGPGLVAPSQSLSITCTVSGFSLTGYGVNWVRQPPGKGLEWLGMIWGDGNTDYNSALKSRLSISKDNSKSQVFLKMNSLHTDDTARYYCARERDYRLDYWGQGTTLTVSS",
        'antigen': "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
    },
    '3HFM': {
        'light':   "DIVLTQSPATLSVTPGNSVSLSCRASQSIGNNLHWYQQKSHESPRLLIKYASQSISGIPSRFSGSGSGTDFTLSINSVETEDFGMYFCQQSNSWPYTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC",
        'heavy':   "DVQLQESGPSLVKPSQTLSLTCSVTGDSITSDYWSWIRKFPGNRLEYMGYVSYSGSTYYNPSLKSRISITRDTSKNQYYLDLNSVTTEDTATYYCANWDGDYWGQGTLVTVSAAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSPRPSETVTCNVAHPASSTKVDKKIVPRDC",
        'antigen': "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
    },
    '5GRJ': {
        'light':   "QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTRVFGTGTKVTVL",
        'heavy':   "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYIMMWVRQAPGKGLEWVSSIYPSGGITFYADTVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARIKLGTVTTVDYWGQGTLVTVSS",
        'antigen': "AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNAALQITDVKLQDAGVYRCMISYGGADYKRITVKVNAPYNKINQRILVVDPVTSEHELTCQAEGYPKAEVIWTSSDHQVLSGKTTTTNSKREEKLFNVTSTLRINTTTNEIFYCTFRRLDPEENHTAELVIPELPLAHPPNER"
    },
    '5Y9J': {
        'light':   "SSELTQDPAVSVALGQTVRVTCQGDSLRSYYASWYQQKPGQAPVLVIYGKNNRPSGIPDRFSGSSSGNTASLTITGAQAEDEADYYCSSRDSSGNHWVFGGGTELTVLGQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS",
        'heavy':   "QVQLQQSGAEVKKPGSSVRVSCKASGGTFNNNAINWVRQAPGQGLEWMGGIIPMFGTAKYSQNFQGRVAITADESTGTASMELSSLRSEDTAVYYCARSRDLLLFPHHALSPWGRGTMVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHHHHHH",
        'antigen': "AVQGPEETVTQDCLQLIADSETPTIQKGSYTFVPWLLSFKRGSALEEKENKILVKETGYFFIYGQVLYTDKTYAMGHLIQRKKVHVFGDELSLVTLFRCIQNMPETLPNNSCYSAGIAKLEEGDELQLAIPRENAQISLDGDVTFFGALKLL"
    },
    '4ETQ': {
        'light':   "QIVLTQSPAIMSAFPGESVTMTCSASSSVSYMYWYQQKPGSSPRLLIYDTSNLASGVPVRFSGSGSGTSYSLTINRLEAEDGATYYCQQWTSYPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNE",
        'heavy':   "QVQLQQSGPELVKPGASVKISCKASGYSFNFYWMHWVKQRPGQGLEWIGMIDPSESESRLNQKFKDKATLTVDRSSSTAHMQLSSPTSEDSAVYYCTRSNYRYDYFDVWGAGTTVTVSSAKTTAPSVYPLAPVCGDTTGSSVTLGCLVKGYFPEPVTLTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVTSSTWPSQSITCNVAHPASSTKVDKKIEPRGG",
        'antigen': "MPQQLSPINIETKKAISNARLKPLDIHYNESKPTTIQNTGKLVRINFKGGYISGGFLPNEYVLSSLHIYWGKEDDYGSNHLIDVYKYSGEINLVHWNKKKYSSYEEAKKHDDGLIIISIFLQVLDHKNVYFQKIVNQLDSIRSANTSAPFDSVFYLDNLLPSKLDYFTYLGTTINHSADAVWIIFPTPINIHSDQLSKFRTLLSLSNHEGKPHYITENYRNPYKLNDDTEVYYSGEIIRAATTSPARENYFMRWLSDLRETLEHHHHHH"
    }
}


# =============================================================================
# Main
# =============================================================================

def main():
    print(f"Device: {DEVICE}")

    # ── load model ────────────────────────────────────────────────────────────
    print("\nLoading TwoStreamSymmetric …")
    model, pkd_bounds = load_model(MODEL_PATH)
    print(f"  pkd_bounds: {pkd_bounds}")

    # Wrap for per-residue IG
    wrapper = TwoStreamPerResidue(model, pkd_bounds).to(DEVICE)
    wrapper.eval()

    # ── load ESM-2 once ───────────────────────────────────────────────────────
    embedder = ESM2PerResidueEmbedder(ESM2_MODEL, DEVICE)

    # ── process each complex ──────────────────────────────────────────────────
    for name, info in COMPLEXES.items():
        print(f"\n{'='*60}")
        print(f"Complex: {name}")
        print(f"{'='*60}")

        print("  Embedding sequences …")
        heavy_emb   = embedder.embed(info['heavy'])    # (L_h, 1280)
        light_emb   = embedder.embed(info['light'])    # (L_l, 1280)
        antigen_emb = embedder.embed(info['antigen'])  # (L_a, 1280)
        print(f"  Heavy: {heavy_emb.shape[0]} residues | "
              f"Light: {light_emb.shape[0]} residues | "
              f"Antigen: {antigen_emb.shape[0]} residues")

        pred = predict_affinity(model, heavy_emb, light_emb, antigen_emb, pkd_bounds)
        print(f"  Predicted pKD: {pred:.3f}")

        # ── Integrated Gradients ──────────────────────────────────────────────
        print(f"  Running IG ({IG_STEPS} steps) …")
        heavy_attr, light_attr, antigen_attr = compute_ig_attributions(
            wrapper, heavy_emb, light_emb, antigen_emb
        )
        print(f"  Attribution ranges: "
              f"Heavy [{heavy_attr.min():.3f}, {heavy_attr.max():.3f}]  "
              f"Light [{light_attr.min():.3f}, {light_attr.max():.3f}]  "
              f"Antigen [{antigen_attr.min():.3f}, {antigen_attr.max():.3f}]")

        # ── figures ───────────────────────────────────────────────────────────
        print("  Generating figures …")

        # 1. Bar plots
        plot_bar(heavy_attr,   info['heavy'],   'Heavy',   name)
        plot_bar(light_attr,   info['light'],   'Light',   name)
        plot_bar(antigen_attr, info['antigen'], 'Antigen', name)

        # 2. Per-chain rectangular heatmaps
        plot_rectangular_heatmap(heavy_attr,   info['heavy'],   'Heavy',   name)
        plot_rectangular_heatmap(light_attr,   info['light'],   'Light',   name)
        plot_rectangular_heatmap(antigen_attr, info['antigen'], 'Antigen', name)

        # 3. Combined 3-panel figure  (Heavy | Light | Antigen)
        plot_combined_heatmaps(
            [heavy_attr, light_attr, antigen_attr],
            [info['heavy'], info['light'], info['antigen']],
            ['Heavy chain', 'Light chain', 'Antigen'],
            name,
            pred_pkd=pred
        )
        print(f"  Saved figures for {name} → {OUTPUT_DIR}/")

    embedder.cleanup()
    print(f"\nDone. All outputs in: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
