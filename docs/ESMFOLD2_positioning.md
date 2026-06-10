# Positioning vs ESMFold2 (and other structure-based predictors)

*How the recent ESMFold2 release relates to our model, how we use it for a fair comparison, and why
our sequence-only affinity predictor remains the right tool for the affinity-regression task.*

## 1. The two models solve different problems

**ESMFold2** (Biohub / EvolutionaryScale, preprint, May 2026) is a **structure-prediction** model: it
turns sequences into atomic-resolution 3D structures of proteins and their complexes, built on ESMC
representations. It is now state of the art for protein–protein and, notably, **antibody–antigen
complex/pose prediction, surpassing AlphaFold3**, and it can run in single-sequence mode for an
order-of-magnitude folding speed-up.

**Our model** solves a different problem: it is a **calibrated binding-affinity regressor** that maps
the three sequences (heavy, light, antigen) directly to a pK_d value. It never predicts a structure.

This distinction matters because **a structure predictor does not output an affinity**. ESMFold2 returns
a 3D pose and per-interface confidence scores (pTM / ipTM), which are measures of *how confident the
model is in the predicted geometry*, not of *how tightly the partners bind*. Confidence and ΔG are only
loosely correlated, so ESMFold2 cannot be read directly as an affinity predictor — it must be coupled to
a downstream scoring step.

## 2. Two principled ways to use ESMFold2 in this study

**(a) As a structure-based affinity baseline (for the comparison table).**
Run ESMFold2 on each complex to obtain a predicted bound structure, then derive an affinity proxy and
correlate it with the measured pK_d under the *same folds* as our model:
- interface confidence (ipTM / interface-pTM),
- a physics/geometry interface score (buried surface area, number of inter-chain contacts/H-bonds, a
  knowledge-based potential), or
- an inverse-folding pseudo-likelihood of the antibody given the predicted interface (e.g. AntiFold /
  ESM-IF), which has been shown to rank affinity-changing mutations zero-shot.

This yields an honest "structure → score" baseline on the SAaIntDB random and cold splits, evaluated by
Pearson/Spearman/RMSE exactly as for our model.

**(b) As a backbone upgrade (ablation / future-proofing).**
Our model consumes frozen **ESM-2 650M** residue embeddings. ESMC (the encoder behind ESMFold2) is a
newer, stronger sequence encoder. Because our interaction module is encoder-agnostic and trains on
cached embeddings, swapping ESM-2 650M → ESMC is a **drop-in change**: re-cache embeddings, retrain the
2.7M-parameter head. We report this as a controlled PLM ablation (alongside ProtBERT, AntiBERTy,
ProGen2), demonstrating the architecture benefits from — but does not depend on — the latest encoder.

## 3. Why our model is the right tool for affinity (the argument to make in the paper)

1. **Direct, calibrated affinity.** We regress pK_d end-to-end in a bounded cosine space; a structure
   pipeline produces a pose plus an *uncalibrated* confidence that must be post-hoc mapped to affinity.
   On the same data we expect our pK_d correlations to **match or exceed** the best structure-derived
   proxy, while being directly interpretable in energy units.
2. **Orders-of-magnitude cheaper.** Inference is a frozen-embedding lookup + a small head — milliseconds
   per complex, CPU-friendly, trivially batched for screening (our cascade-screening experiments rely on
   this). Even in single-sequence mode, ESMFold2 runs a 3D diffusion model per complex and then needs a
   separate scoring pass; for ranking millions of designs this gap is decisive.
3. **No MSA, no structure required.** We need only the three sequences, so the method applies to
   *de novo* designs, nanobodies (n = 662 here), and deep-mutational landscapes (zero/few-shot) where no
   co-structure exists. ESMFold2 is strongest with evolutionary information (MSA); our model deliberately
   avoids that dependency.
4. **Robust on the regimes that matter.** We report cold (no-PDB-overlap), cold-antigen and
   30%-sequence-identity splits, and nanobody/antibody stratification — the out-of-distribution settings
   relevant to prospective discovery, where confidence-based proxies are least reliable.

**Framing sentence for the paper.** *"ESMFold2 sets a new bar for antibody–antigen structure prediction,
but structure confidence is not affinity. We therefore include a structure-based affinity baseline
(ESMFold2 pose → interface score) and show that our sequence-only regressor attains equal or better pK_d
correlation at a fraction of the computational cost, while extending to nanobodies and low-data
mutational landscapes where no co-structure is available."*

## 4. Minimal experiment to substantiate the comparison

1. Predict complexes for the SAaIntDB benchmark set (and the 5 explainability complexes) with ESMFold2
   (single-sequence mode; MSA mode as an upper bound).
2. Score each with (i) ipTM, (ii) interface contacts/BSA, (iii) AntiFold pseudo-likelihood.
3. Correlate each proxy with measured pK_d on the **same random and cold folds**; report
   Pearson/Spearman/RMSE + wall-clock per complex.
4. Add one row per proxy to the SOTA comparison table next to "Ours" and the sequence baselines
   (MVSF-AB, etc.), plus a cost column (params, s/complex, MSA yes/no).

This keeps the comparison apples-to-apples and lets the figure make the cost/accuracy trade-off explicit.

## Sources
- [Biohub — a world model of protein biology (ESMFold2 release)](https://biohub.org/news/world-model-of-protein-biology/)
- [Nature news — open-source model predicts 1 billion protein structures](https://www.nature.com/articles/d41586-026-01686-3)
- [Latent Space — ESMFold2: the Bitter Lesson for proteins (Alex Rives, Biohub)](https://www.latent.space/p/esmfold2)
- [AntiFold: improved structure-based antibody design using inverse folding (Bioinformatics Advances)](https://academic.oup.com/bioinformaticsadvances/article/5/1/vbae202/8090019)
