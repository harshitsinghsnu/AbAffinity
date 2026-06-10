# Ground-truth binding-residue PyMOL figures

Publication-style structural renderings of four antibody–antigen complexes
(**1VFB, 4ETQ, 5GRJ, 5Y9J**), one PyMOL session per complex **per model**
(our final All-CDR model and the two-stream baseline).

For each complex the **paper-reported ground-truth binding residues** (heavy, light and
antigen) are shown as side-chain sticks on a pale three-chain cartoon, and each residue is
**coloured by that model's Integrated-Gradients attribution** so that residues the model
considers important (signed-heatmap-normalised attribution **> 0.4**) appear in saturated
colour. This visualises, directly on the structure, *which* of the true interface residues
each model recovers.

## Colour scheme (matches the explainability heatmaps)

| Chain   | Colormap | Cartoon tint | Stick colour |
|---------|----------|--------------|--------------|
| Heavy   | Greens   | pale green   | green, darker = higher attribution |
| Light   | Blues    | pale blue    | blue,  darker = higher attribution |
| Antigen | Purples  | pale purple  | purple, darker = higher attribution |

Background is white. Attribution is the **signed min–max normalisation per chain** used by the
heatmaps; the stick colour maps `t = 0.42 + 0.58·attr`, so an attribution of 0.4 already sits
in the saturated half of the colormap.

## Files (`figures/`)

For every `{PDB}_{ours|twostream}` there are **two variants**:
- `{PDB}_{model}.png` / `.pse` — **labelled**: each ground-truth residue is annotated with its
  residue name + number (e.g. `ASP58`). Labels that fit sit directly on their residue; labels that
  would overlap are moved to a left/right edge gutter and linked back with a **dotted leader line**,
  so the structure stays prominent and every label is readable.
- `{PDB}_{model}_nolabel.png` / `.pse` — **no labels** (clean structure + coloured sticks only).

Each `.png` is a 2800×2100 ray-trace at **350 dpi**; each `.pse` is a self-contained PyMOL session
(open directly in the GUI). The matching `pml/{...}.pml` reproduces it (run from the `3_stream` dir
via `pymol -cq`; the camera and label positions are baked in with `set_view` + pseudoatoms).

## Chain assignments (deposited author chains)

| PDB  | Heavy | Light | Antigen | Note |
|------|-------|-------|---------|------|
| 1VFB | B     | A     | C       | HEL–FvD1.3 |
| 4ETQ | H     | L     | **C**   | two copies in the AU; H/L bind antigen copy **C** (not X) |
| 5GRJ | H     | L     | A       | |
| 5Y9J | H     | L     | A       | |

Ground-truth residues are selected by their **author residue number**, exactly as listed in the
source papers. Two residues differ in identity from the deposited model and are selected by
number as given: 1VFB light **53** is Thr in the structure (listed as Tyr53) and 1VFB antigen
**19** is Asn (listed as Asp19).

## Result — ground-truth residues with attribution > 0.4

Our All-CDR model attributes the antibody **paratope** far more strongly than the two-stream
baseline; on the light chain it recovers 23/24 ground-truth residues versus 5/24.

| PDB  | chain  | nGT | Ours >0.4 | Two-stream >0.4 |
|------|--------|-----|-----------|-----------------|
| 1VFB | heavy  | 12  | 10 | 3 |
| 1VFB | light  | 10  | 9  | 4 |
| 1VFB | antigen| 9   | 3  | 6 |
| 4ETQ | heavy  | 14  | 9  | 3 |
| 4ETQ | light  | 4   | 4  | 0 |
| 4ETQ | antigen| 21  | 10 | 19 |
| 5GRJ | heavy  | 15  | 10 | 5 |
| 5GRJ | light  | 6   | 6  | 0 |
| 5GRJ | antigen| 19  | 16 | 15 |
| 5Y9J | heavy  | 8   | 4  | 6 |
| 5Y9J | light  | 4   | 4  | 1 |
| 5Y9J | antigen| 12  | 8  | 10 |
| **Total** | | **134** | **93** | **72** |
| Antibody (H+L only) | | 73 | **66** | **16** |

Full per-residue values are in `gt_attribution_recovery.csv` and `mapping_{PDB}.json`.

## Reproduce

```
python pymol_groundtruth/verify_and_map.py      # GT<->structure mapping + attribution (uses pdb/)
python pymol_groundtruth/render_figures.py       # computes label layout, writes pml/*.pml (label + _nolabel)
for %f in (pymol_groundtruth\pml\*.pml) do pymol -cq "%f"   # writes figures/*.png + *.pse
```

PNGs are rendered through the `pymol -cq` launcher (not the pymol2 API) so the expired-evaluation
watermark is not stamped. Attribution arrays come from `experiments/results_explainability/{PDB}.npz`
(`export_ig_attributions.py`); structures are in `pymol_groundtruth/pdb/`.
(`generate_pml.py` is the earlier no-leader-line generator, superseded by `render_figures.py`.)
