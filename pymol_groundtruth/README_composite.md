# Composite interface figures (reference-image style)

Publication composites in the style of the supplied reference: a central **antibody (salmon) +
antigen (green)** overview with transparent surface + cartoon, **dashed circles** on the interface
sub-regions, and **3 zoom-in inset boxes** (connected by black leader lines) showing the ground-truth
interface residues as **sticks**, **hydrogen bonds as yellow dashes**, and **RES-NUM labels**.

- Residue sticks are shaded by **that model's IG attribution** within the partner colour
  (antibody = Reds ramp, antigen = Greens ramp), so `ours` and `twostream` differ.
- Interface sub-regions are found by k-means (k=3) on the ground-truth residue Cα coordinates;
  each cluster becomes one inset.
- Hydrogen bonds are N/O atom pairs < 3.5 Å between distinct ground-truth residues (yellow dashes).
- Labels inside an inset are placed on their residue when possible, else nudged to the inset edge
  with a grey dotted leader (same declutter as the single-view figures).

## Files (`figures_composite/`)
- `{PDB}_{ours|twostream}_composite.png` / `.pdf` — the assembled composite (300 dpi).
- intermediate renders `{PDB}_{model}_overview.png` and `_inset0..2.png` (PyMOL ray traces).

PNGs are rendered through `pymol -cq` (no evaluation watermark); the matplotlib step assembles the
overview + insets, draws the dashed circles and connector lines.

## Reproduce
```
python pymol_groundtruth/composite_figures.py
```
Uses `mapping_{PDB}.json` (GT residues + attribution), structures in `pymol_groundtruth/pdb/`,
writes `.pml` to `pymol_groundtruth/pml_composite/`. Chain assignments and the 1VFB Thr53/Asn19
caveats are the same as in the main [README.md](README.md).
