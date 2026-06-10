# Methods

*Target journal: Nature Communications. Written as continuous Methods prose — no source-code or
function references, and only the equations that are central to the model are made explicit; all other
operations are described conceptually.*

## Datasets and affinity representation

The primary dataset is **SAaIntDB**, a curated collection of **2,575 structurally resolved
antibody–antigen complexes**, each comprising an antibody heavy chain, an antibody light chain (absent
for single-domain antibodies), and an antigen, together with an experimentally measured equilibrium
dissociation constant and antigen-identity annotations. Single-domain antibodies (nanobodies/VHH;
n = 662) are retained throughout: because they lack a light chain, the heavy-chain representation is
also placed in the light-chain position so that the three-stream interface is defined consistently for
every complex. This design lets a single model span conventional paired antibodies and nanobodies
without architectural changes.

All affinities are expressed on the pK_d scale,

  pK_d = −log₁₀(K_d),

so that larger values denote stronger binding. Reporting in pK_d gives a bounded, monotonic target
that is well matched to a similarity-based regression objective and is consistent with recent
sequence-based affinity benchmarks.

For comparison with prior sequence-based predictors we additionally evaluate on four widely used
collections derived from the same public release used by earlier work: a natural-complex set (SAbDab,
578 complexes), two mutational sets (AB-Bind, 1,089 variants; SKEMPI 2.0, 387 variants), and an
independent held-out benchmark of 38 natural complexes. Free-energy values for these sets are
converted to pK_d through the standard thermodynamic relation ΔG = RT ln K_d at 298.15 K. To probe
low-data adaptation we use seven single-antigen deep-mutational landscapes, evaluated in zero-shot and
few-shot regimes.

## Protein-language-model features and CDR-focused pooling

Sequences are encoded with a frozen **ESM-2 (650M-parameter) protein language model**; its weights are
never updated, so the entire learning problem reduces to training a small interaction module
(~2.7 million parameters) on top of fixed, information-rich residue embeddings. Working from cached
embeddings makes large-scale cross-validation and ablation tractable and guarantees that every
architectural comparison consumes identical upstream features.

A central methodological choice is **how each chain is summarised into a fixed-length vector**. The
antigen and the light chain are reduced by averaging their residue embeddings, which preserves global
sequence context. For the **heavy chain**, however, specificity is concentrated in three short
complementarity-determining regions (CDRs), and averaging over the entire chain dilutes this signal
with conserved framework residues. We therefore introduce **All-CDR pooling**: the heavy chain is
numbered under the IMGT convention to locate CDR-H1, CDR-H2 and CDR-H3, and only the residue embeddings
within these loops are averaged to form the heavy-chain representation. When the numbering procedure
fails, a motif-based fallback identifies CDR-H3, and if no loop can be resolved the chain reverts to
full-length averaging. Conceptually, All-CDR pooling acts as a biologically grounded attention prior:
it concentrates the antibody representation on the paratope loops that determine binding while leaving
the antigen description global, because the relevant epitope is not known a priori. This focusing
improves accuracy on the primary dataset and transfers better to unseen mutational landscapes than
full-chain pooling.

Each complex therefore enters the network as three fixed-length vectors — heavy, light and antigen.

## Model architecture

The model is a **three-stream factorised network** that keeps the heavy chain, the light chain and the
antigen as separate but interacting representations until the final affinity readout. The design is
motivated by antibody biology: heavy and light chains arise from independent recombination events and
contribute unequally to binding, and the antigen should be treated as an active partner whose
binding-relevant features must be distinguished from generic protein context.

**Stream-specific projection.** Each raw embedding is first normalised and then mapped, through its own
learnable projection, into a shared 256-dimensional interaction space. Separate projections let the
model learn chain-specific adapters before any cross-molecular comparison, and the initial
normalisation prevents any single stream from dominating purely through differences in embedding scale.

**Heavy–light coupling and adaptive fusion.** The projected heavy and light vectors are allowed to
exchange information through self-attention, capturing the cooperativity of the two chains within the
paratope. They are then combined into a single antibody representation by a **learned fusion gate**
that produces two normalised weights and forms a convex combination,

  a_fused = w_H·h + w_L·l,   (w_H, w_L) = softmax(·),

so that the model can adaptively emphasise the heavy or the light chain on a per-complex basis. This is
biologically meaningful because some interfaces are heavy-chain dominated whereas others rely more on
the light chain, and for nanobodies the gate naturally collapses onto the single available chain.

**Gated cross-attention as a denoising interaction.** Antibody and antigen representations interact
through a **gated cross-attention** block. Rather than computing a token-by-token attention matrix, the
block forms a *feature-wise multiplicative* match between the antibody query and the antigen
key/value, and modulates it with a learned, query-dependent gate. With learnable projections of the
query and key/value, the core operation is

  interaction = tanh(W_q q ⊙ W_k k),   g = σ(W_g q),
  output = LayerNorm( q + W_o[ (interaction ⊙ W_v k) ⊙ g ] ),

where ⊙ denotes element-wise multiplication and σ the logistic function. The tanh interaction
emphasises latent dimensions in which antibody and antigen features are jointly informative, while the
sigmoid gate **selectively preserves or suppresses individual dimensions** according to the binding
context. Conceptually the gate acts as a context-dependent denoiser that filters out components of the
language-model embedding reflecting fold or family identity rather than binding. The block is applied
bidirectionally and stacked (two interaction layers), with the fused antibody conditioning the antigen
along the principal prediction pathway. Disabling the gate reduces the module to an ungated bilinear
interaction and, as shown by our ablations, sharply degrades generalisation — establishing that the
gate is functionally essential rather than cosmetic.

**Similarity-based readout.** The interacting antibody and antigen representations are passed through
separate output heads into partner-specific subspaces, L2-normalised, and compared by cosine
similarity; the cosine score is the predicted affinity. Predicting in a bounded similarity space —
rather than through an unconstrained regression head — keeps the two partners geometrically comparable,
reduces representational collapse and, empirically, transfers better to unseen complexes than MLP or
tree-based heads on the same features.

**Baselines.** Two code-matched baselines isolate the contribution of explicit chain factorisation. A
**two-stream** model fuses the heavy and light chains into a single antibody embedding (joining the two
sequences with explicit boundary tokens before encoding) but retains the same gated interaction and
cosine readout, so the only missing ingredient is heavy/light disentanglement. A **concatenation +
multilayer-perceptron** baseline regresses affinity directly from the concatenated heavy, light and
antigen embeddings, with no partner conditioning or similarity-based objective.

## Training objective and optimisation

During training the experimental pK_d values are linearly mapped to the interval [−1, 1] using bounds
computed on the training split, and the model is optimised to minimise the mean squared error between
this rescaled target and the predicted cosine similarity. This couples the regression target directly
to the geometry of the shared embedding space. Optimisation uses AdamW (learning rate 1×10⁻⁴, weight
decay 0.01) with gradient clipping and early stopping on validation correlation; the language-model
backbone remains frozen, and only the projection, fusion, interaction and head parameters are learned.

## Evaluation and statistical analysis

Performance is assessed by **ten-fold cross-validation** under two complementary splitting protocols
that probe different notions of generalisation.

In the **random split**, the complete set of complexes is partitioned into ten equally sized folds by
shuffling and dividing at the level of individual complexes; each fold is held out once for testing
while the remaining nine are used for training. Because the partition is made over individual
complexes, two complexes that originate from the *same* deposited structure (the same PDB entry) — for
example several antibody variants co-crystallised with a common antigen — may fall on opposite sides of
the train/test boundary. The random split therefore measures interpolation performance within the
distribution of the dataset and is the standard setting for comparing architectures and pooling
strategies.

In the **cold (structure-grouped) split**, the partition is instead made over **unique PDB
identifiers** rather than over individual complexes: the set of distinct structures is shuffled and
divided into ten groups, and every complex belonging to a held-out structure is assigned *entirely* to
the test fold. This guarantees that **no PDB entry is shared between training and test** in any fold, so
that all complexes derived from a given structure (and hence closely related antibody–antigen pairs)
are evaluated only as unseen examples. The cold split is substantially more demanding than the random
split because the model cannot exploit near-duplicate complexes seen during training, and it directly
estimates generalisation to structurally novel targets — the regime most relevant to prospective
antibody discovery.

For the **external mutational datasets (AB-Bind and SKEMPI)** we additionally evaluate two stricter
hold-out protocols. In the **cold-antigen** protocol, the antigen PDB code is parsed from each complex
identifier and a fraction of the distinct antigen structures is held out entirely for testing, so that
no antigen PDB is shared between training and test. In the **sequence-similarity** protocol, antibody
sequences (heavy and light concatenated) and antigen sequences are independently clustered by greedy
single-linkage agglomeration at a 30% pairwise sequence-identity threshold; a fraction of the antibody
clusters and of the antigen clusters is held out, and any complex whose antibody *or* antigen belongs
to a held-out cluster is placed in the test set. This removes both near-identical structures and
near-identical sequences from the training distribution, providing the most conservative estimate of
out-of-distribution generalisation. Both protocols are run over three seeds.

Every experiment is repeated over three random seeds, which control both the fold assignment and the
network initialisation, and results are reported as mean ± standard deviation across seeds using
Pearson correlation, Spearman rank correlation and root mean squared error (in pK_d).

Uncertainty on the headline predictions is quantified by a **non-parametric bootstrap** (1,000
resamples of the held-out predictions, reported as 95% percentile intervals) and a **label-permutation
test** (1,000 permutations) for statistical significance. To distinguish global ranking from
within-target ranking we additionally compute **per-antigen** correlations for targets with at least
three complexes and aggregate them in Fisher-z space, reporting the average within-target correlation
alongside the pooled value; this separates the model's ability to rank diverse complexes from its
ability to rank antibodies against a common antigen. Finally, performance is stratified into
conventional paired antibodies and single-domain nanobodies to confirm that both regimes are handled
comparably.

## Ablation studies

A series of controlled ablations isolates the contribution of each design choice: (i) **pooling** —
full-chain averaging versus All-CDR pooling (heavy-only and heavy + light); (ii) **architecture** —
three-stream versus the fused two-stream and concatenation baselines; (iii) **language-model backbone**
— the same architecture paired with alternative encoders, including a mixed antibody-specific/general
encoder configuration; (iv) **gating** — evaluating the trained model with the gate enabled, fully
open, fully closed, fixed, or randomised, to quantify its effect while holding all other weights fixed;
and (v) **antigen dependence** — replacing the antigen with zeroed, averaged or shuffled embeddings to
test whether predictions genuinely use partner information rather than antibody-only regularities.

## Interpretability by integrated gradients

To test whether the model attends to physically meaningful determinants, we compute residue-level
attributions with **integrated gradients**. Because the network consumes pooled representations, it is
wrapped in a differentiable pooling layer so that gradients propagate back to individual residues. For
each residue the attribution integrates the model's sensitivity along a straight path from a zero
baseline to the observed embedding,

  A(x_i) = (x_i − x̃_i) · ∫₀¹ ∂f(x̃ + α(x − x̃))/∂x_i dα,

evaluated with fifty integration steps and summed over the embedding dimension to give one signed score
per residue, then min–max normalised within each chain for visualisation. Attributions are computed for
the final All-CDR model and, for comparison, the two-stream baseline, on five structurally
characterised complexes, and rendered as per-chain residue-importance maps with a fixed colour
assignment for the heavy chain, light chain and antigen. We summarise localisation by a CDR-enrichment
score — the fraction of total heavy-chain attribution falling within the CDR loops, divided by the
fraction of heavy-chain residues that are CDRs — where values above one indicate preferential
attribution to the paratope.

## Reporting summary

ESM-2 embeddings are precomputed once and cached so that all comparisons use identical upstream
features. Software versions are pinned in the accompanying environment specification (PyTorch
2.2.1 + CUDA 12.1, Transformers 4.47.1, captum 0.8.0, ANARCI 2026.2, scikit-learn 1.6.1, NumPy 1.26,
pandas 2.3, Matplotlib 3.10) to ensure reproducibility.
