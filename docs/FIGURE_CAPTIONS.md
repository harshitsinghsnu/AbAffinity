# Figure captions (publication notebook)

Captions for the figures assembled in `notebooks/BALM_AbAg_Panels.ipynb`. All panels are 400 DPI
(PDF + PNG); the model is the final **All-CDR three-stream** model; the dataset is **SAaIntDB**
(2,575 antibody–antigen complexes, including 662 nanobodies) unless otherwise stated. Correlations are
Pearson *r* and Spearman ρ; error is in pK_d.

**Figure 1 | Affinity distribution of the benchmark datasets.** Distribution of experimental binding
affinity (pK_d) across SAaIntDB and the external datasets, showing the dynamic range and class balance
spanned by the evaluation. (`pubfig1_pkd_distributions`)

**Figure 2 | Three-stream All-CDR model outperforms fused and concatenation baselines.** Ten-fold
cross-validation on SAaIntDB for the proposed three-stream All-CDR model, a fused two-stream baseline,
and a concatenation+MLP baseline, under both the random and the cold (no-PDB-overlap) split, for
Pearson *r*, Spearman ρ and RMSE (six panels). Explicit heavy/light factorisation gives the best
correlation and lowest error in every setting; significance between models is annotated.
(`pubfig2_architecture_saaintdb`)

**Figure 3 | Robustness of the All-CDR model.** Bootstrap-resampled performance (1,000 resamples) for
Pearson, Spearman and RMSE, with per-PDB and per-antigen distributions of correlation, demonstrating
that performance is stable across resamples and across individual targets. (`pubfig3_robustness`)

**Figure 3b | Overall versus target-level (Fisher-aggregated) ranking.** Per-target correlations
aggregated in Fisher-z space compared with the pooled (overall) correlation, for both per-antigen and
per-PDB groupings (Pearson and Spearman). The model ranks complexes well overall and retains positive
within-target ranking, with the expected larger variance across individual targets.
(`pubfig3b_overall_vs_fisher`)

**Figure 3c | Overall versus intra-target performance across all metrics.** Direct comparison of
overall (pooled) versus intra-target (per-antigen, Fisher-aggregated for correlations) Pearson,
Spearman and RMSE; whiskers denote ±1 s.d. across antigen targets. (`pubfig3c_overall_vs_intratarget`)

**Figure 4 | ESM-2 is the strongest sequence encoder.** Ten-fold CV on SAaIntDB with the same
architecture paired with ESM-2, ProGen2, AntiBERTy and ProtBERT, plus a mixed antibody-specific/general
configuration (AntiBERTy heavy/light + ESM-2 antigen), for all three metrics. ESM-2 gives the best
overall transfer. (`pubfig4_plm_comparison`)

**Figure 5 | Pooling ablation.** Full-chain mean pooling versus All-CDR (heavy + light) versus All-CDR
(heavy-only) on SAaIntDB (random, three seeds), for all three metrics. CDR-focused pooling improves
over full-chain pooling; heavy-only and heavy+light are statistically indistinguishable, and heavy-only
is adopted as the final model for simplicity and nanobody compatibility. (`pubfig5_pooling_ablation`)

**Figure 5b | Headline pooling comparison.** Focused comparison of the final All-CDR (heavy
H1+H2+H3) model against full-chain mean pooling across Pearson, Spearman and RMSE. (`pubfig5b_meanpool_vs_allcdr`)

**Figure 6 | Zero-shot and few-shot transfer to mutational landscapes.** Zero-shot prediction and
few-shot fine-tuning with 10%, 20% and 30% labelled data (three seeds, error bars) on the
single-antigen mutational datasets, for Pearson and Spearman, with a summary table at 30%. Modest
labelled data substantially improves ranking over zero-shot. (`pubfig6_transfer`)

**Figure 7 | Top-K recovery.** Fraction of the true top-ranked binders recovered within the model's
top-K predictions, quantifying prioritisation performance for screening. (`pubfig7_topk`)

**Figure 8 | Saturation mutagenesis on a CDR-H3 landscape.** Predicted mutational-effect heatmap over a
heavy-chain CDR-H3 saturation scan, with Pearson, Spearman and RMSE versus measured effects in the
panel titles; charged-residue substitutions produce the largest predicted shifts. (`pubfig8_mutagenesis`)

**Figure 9 | Gating ablation (all metrics).** Post-hoc gate ablation on SAaIntDB with the trained model
evaluated using the learned gate, a fixed gate (0.5), a fully open gate (g=1), a closed gate (g=0) and
a randomised gate, for Pearson, Spearman and RMSE. Disabling the gate collapses performance
(learned *r* = 0.94 → closed *r* = 0.61). (`pubfig9_gating_saaintdb`)

**Figure 10 | Amino-acid substitution impact.** Mean true ΔpK_d, mean predicted ΔpK_d, and prediction
error across the 20×20 wild-type→mutant substitution matrix (residues coloured by physicochemical
class), showing close agreement on the dominant substitution trends. (`pubfig10_aa_mutation_impact`)

**Figure 11 | Comparison with MVSF-AB.** Dataset-wise comparison of the All-CDR model against the
MVSF-AB baseline on the natural and mutational datasets (Pearson and RMSE; MVSF-AB reported only these),
with seed-level error bars on our model. Our model exceeds MVSF-AB on every setting. (`pubfig11_vs_mvsf`)

**Figure 12 | Predicted versus measured affinity (SAaIntDB).** Out-of-fold predicted versus measured
pK_d for the All-CDR model on SAaIntDB, with Pearson, Spearman and RMSE annotated. (`pubfig12_scatter`)

**Figure 13 | Performance on antibodies and nanobodies.** Overall versus paired-chain antibody versus
single-domain nanobody (VHH) performance (Pearson, Spearman, RMSE), showing comparable accuracy on both
antibody classes. (`pubfig13_antibody_vs_nanobody`)

**Figure 14 | Per-antigen intra-target ranking.** Per-antigen Pearson and Spearman correlation for the
top targets (n ≥ 3 complexes), demonstrating the model's ability to rank antibodies against a common
antigen. (`pubfig14_per_antigen_ranking`)

**Figure 15 | Bootstrap confidence and significance.** Bootstrap distributions (1,000 resamples) for
Pearson, Spearman and RMSE with shaded 95% confidence bands; a label-permutation test gives p < 0.001.
(`pubfig15_bootstrap`)

**Figure 17 | Predicted versus measured affinity with 95% confidence band (SAaIntDB).** Out-of-fold
scatter with an ordinary-least-squares fit and a bootstrap 95% confidence band for the regression line.
(`pubfig17_scatter_ci_band`)

**Figure 18 | Predicted versus measured affinity on external datasets.** Out-of-fold scatter of the
All-CDR model on SAbDab, AB-Bind and SKEMPI, each with an OLS fit and 95% confidence band and metrics
annotated. (`pubfig18_scatter_natural`)

**Figure 19 | The gating mechanism and its impact.** **a** Learned antibody→antigen gate σ(·) averaged
over all SAaIntDB complexes, reshaped to a 16×16 latent grid, revealing a sparse, structured pattern
that amplifies a subset of dimensions and suppresses others. **b** Distribution of per-dimension gate
values (≈61% amplified, ≈6% suppressed). **c** Post-hoc gate ablation: disabling the gate reduces
Pearson *r* from 0.94 to 0.61 (−0.33), confirming the gate is functionally essential.
(`pubfig19_gating_impact`)

**Supplementary Figure S1 | Stream-intervention and per-CDR-region analyses.** Performance under
zeroing of the light, heavy or antigen stream and under mean/shuffled antigen embeddings (Pearson and
Spearman), together with per-CDR-region (H1/H2/H3) performance across all three metrics, quantifying the
model's dependence on authentic antigen information and on each CDR loop. (`pubfigS1_supplementary`)

---

# Composite panels (`BALM_AbAg_Panels.ipynb`)

Self-contained Nature-style multi-panel figures (400 DPI PDF+PNG). Bar/curve panels embed their exact
numbers; scatter and heatmap panels read raw data. "Ours" uses one fixed colour (vermillion) throughout;
comparators use the colour-blind-safe Okabe-Ito palette; ordered series use sequential colormaps.

**Panel 1 | Datasets and pK_d distributions.** **a** Normalised affinity (pK_d) distributions across
SAaIntDB and the four external datasets. **b** Dataset sizes (log scale). **c** SAaIntDB composition:
paired antibodies (1,913) versus nanobodies (662). **d** pK_d dynamic range (mean marked) per dataset.
(`panel1_data`)

**Panel 2 | Main comparison.** **a** SAaIntDB architecture comparison — Ours (All-CDR), two-stream and
concat+MLP under random (solid) and cold (hatched) splits, for Pearson, Spearman and RMSE. **b** Antibody
versus nanobody performance, one subplot per metric (Pearson, Spearman, RMSE), with error bars computed
across the ten cross-validation folds. **c** Ours versus MVSF-AB across SAbDab, AB-Bind, SKEMPI and the
Benchmark for Pearson and RMSE (MVSF-AB did not report Spearman). **d** Overall (pooled) versus
intra-target performance across all three metrics. *Intra-target* metrics are obtained by computing the
correlation **within each antigen target** (antigens with ≥3 complexes), then averaging the per-target
correlations across targets in Fisher-z space — r̄ = tanh(mean over targets of arctanh r_i) — which
weights each target equally regardless of size; error bars are ±1 s.d. across targets for intra-target
and ±1 s.d. across the ten folds for overall. Intra-target therefore measures the model's ability to
**rank antibodies against a common antigen**, the regime relevant to affinity maturation, whereas
overall measures ranking across the whole heterogeneous pool. (`panel2_main`)

**Panel 3 | Protein-language-model comparison.** Ten-fold CV on SAaIntDB for ESM-2, ProtBERT, AntiBERTy,
ProGen2 and the mixed AntiBERTy(H/L)+ESM-2(Ag) configuration, for Pearson, Spearman and RMSE (sequential
colours). ESM-2 is strongest. (`panel3_plm`)

**Panel 4 | Zero-shot and few-shot transfer.** One subplot per mutational dataset (1mlc, aayl51, 4fqi,
koenig, trastuzumab, warszawski). In each, the dotted line is the zero-shot Pearson, the solid line with
error bars is few-shot Pearson (mean ± s.d. over 3 seeds) and the dashed line is few-shot Spearman, as a
function of labelled-data fraction (10/20/30%). (`panel4_transfer`)

**Panel 5 | Gating ablation, antigen-dependence and prioritisation.** **a** Post-hoc gate ablation
(learned/fixed/open/random/closed) for Pearson, Spearman and RMSE — closing the gate collapses
performance. **b** Stream-intervention controls (zero/mean/shuffled streams) across all three metrics
(Pearson, Spearman, RMSE; error bars are ±1 s.d. across folds). **c** Per-CDR-region (H1/H2/H3)
performance across all three metrics. **d** **Top-K binder recovery** on the AAYL51 mutational test set:
the fraction of the experimentally top-ranked binders that fall within the model's own top-K% of
predictions, for K = 5/10/20/30%, compared against a mean-pool model, a zero-shot model, and a random
baseline. *How it works:* variants are ranked by predicted affinity, the top-K% are taken, and recovery
= (number of true top-K% binders captured) / (size of the top-K% set). *Why it matters:* in real
campaigns only a small fraction of designs can be assayed, so the practically relevant question is not
the global correlation but **how many genuine strong binders survive an aggressive top-K cut** — i.e.
the model's value as a triage/prioritisation engine. Across all cut-offs the All-CDR model recovers the
most true top binders (e.g. **51.7% at the top-30% cut versus 30% for random**), and it exceeds the
zero-shot and mean-pool baselines, demonstrating useful early-discovery prioritisation. (`panel5_gating_supp`)

**Panel 6 | Mutational heatmaps.** **a** Saturation-mutagenesis measured ΔpK_d for CDR-H1/H2/H3, with the
x-axis labelled by wild-type residue and position and amino-acid labels coloured by physicochemical
class. **b** 20×20 wild-type→mutant substitution matrix (true ΔpK_d, predicted ΔpK_d, and error), with
amino-acid labels coloured by class (charged red, polar blue, hydrophobic green, special purple).
(`panel6_heatmaps`)

**Panel 7 | Predicted versus measured affinity.** Two-row layout of scatters (best representative
fold/seed) with OLS fit and bootstrap 95% confidence band for **a** SAaIntDB, **b** SAbDab, **c** AB-Bind,
**d** SKEMPI and **e** the Benchmark set, with Pearson, Spearman, RMSE and n annotated. (`panel7_scatters`)
