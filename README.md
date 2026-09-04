<h1 align="center">bioIOT <span style="font-size:60%">(R)</span></h1>

<p align="center">
  <b>Inverse optimal transport for single-cell state transitions and pseudotime</b>
</p>

<p align="center">
  <a href="https://github.com/XTSgreen/bioIOT-R/actions/workflows/ci.yml"><img src="https://github.com/XTSgreen/bioIOT-R/actions/workflows/ci.yml/badge.svg" alt="R-CMD-check"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="License"></a>
</p>

---

<p align="center"><a href="README.md"><b>English</b></a> | <a href="https://github.com/XTSgreen/bioIOT-R/blob/main/README.zh-CN.md">简体中文</a></p>

**bioIOT** is an R implementation of semi-relaxed inverse optimal transport
(IOT) for single-cell trajectory analysis. Given state-transition features,
source/target state masses and observed transitions, it learns feature
weights θ such that the soft-marginal OT plan induced by the linear cost
`C = -einsum(φ, θ)` reproduces the data — and turns the fit into state
transition matrices, random-walk pseudotime and ggplot2 visualisation.

The solver was developed and validated as part of a research project on
treatment-resistance state transitions, and is packaged here for general
use.

## Why bioIOT?

- **Identifiable by construction.** Hard-marginal OT makes pure column
  features unidentifiable; bioIOT's KL-soft column anchoring restores
  identifiability while preserving the target-composition constraint.
- **Exact implicit gradients.** Anderson-accelerated fixed-point solve in
  the forward pass; the implicit function theorem in the backward pass —
  numerically stable where unrolled backpropagation diverges.
- **Self-contained.** The solver is pure base R; the only runtime dependency
  beyond it is ggplot2. Seurat / SingleCellExperiment are soft-gated.
- **End-to-end.** Cell embedding in, transition matrix + pseudotime + plots
  out.

## Installation

```r
# from GitHub
remotes::install_github("XTSgreen/bioIOT-R")
```

<details>
<summary>From a local clone</summary>

```r
install.packages(".", repos = NULL, type = "source")
# or: R CMD build . && install.packages("bioIOT_0.2.0.tar.gz", repos = NULL)
```

</details>

## Quick start

```r
library(bioIOT)

# 1) A reproducible synthetic dataset with ground truth
sim <- simulate_iot_states(K = 6, seed = 1)   # or data(demo_iot_states)

# 2) Fit feature weights: exact implicit gradients + two-stage debias +
#    multi-restart. Scenarios may differ in the number of states.
fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T_true, n_restart = 2)
fit; summary(fit)

# 3) Trajectory layer
Q  <- transition_matrix(fit)                        # (K, K) transitions
pt <- pseudotime_from_transition(Q, root = "S1")    # random-walk pseudotime

# 4) Visualisation
plot_transition_heatmap(Q)                          # annotated heatmap
plot_transition_flow(Q, sim$embedding)              # CellRank-style arrows
plot_theta(fit)                                     # weights + support

# 5) Straight from cell-level objects
res <- runIOT(sim$cell_embedding, sim$cell_state,
              from = sim$cell_time == "t0", to = sim$cell_time == "t1",
              root = "S1")
# SingleCellExperiment: runIOT(sce, state_col = "state", time_col = "time",
#                              from = "t0", to = "t1", dimred = "PCA")
# Seurat: runIOT(obj, group.by = "state", split.by = "time",
#                 from = "t0", to = "t1", reduction = "pca")
```

Without observed transitions `T_obs`, `runIOT()` solves the plan with
uniform feature weights; with `T_obs` (e.g. from lineage or clone data) it
fits the weights first.

## Showcase: HSMM differentiation (real public data)

The [showcase directory](https://github.com/XTSgreen/bioIOT-R/tree/main/showcase)
runs the full pipeline on the human
skeletal-muscle myoblast time course (HSMMSingleCell, 271 cells × 47k genes,
0/24/48/72 h; Shin et al. 2015) and benchmarks it against Slingshot:

| Pseudotime metric (Spearman) | value |
|---|---|
| bioIOT pseudotime ~ known Hours (cells) | 0.251 |
| Slingshot pseudotime ~ known Hours (cells) | 0.253 |
| bioIOT state pseudotime ~ Slingshot state pseudotime | 0.600 |

With 30 cells of sampling noise per state, `fit_iot` recovers the true
transition matrix ~3× more accurately than the raw noisy transitions
(mean per-row L1 0.013 vs 0.039; 50 replicates). Provenance and licensing of
all showcase data are documented in
[showcase/LICENSE_AUDIT.md](https://github.com/XTSgreen/bioIOT-R/blob/main/showcase/LICENSE_AUDIT.md).

## API overview

| Layer | Functions |
|---|---|
| Core solver | `soft_sinkhorn()`, `row_conditional()`, `make_cost()`, `row_ce_loss()`, `zscore_phi()` |
| Fitting | `fit_iot()` (with `print`/`summary` methods) |
| Trajectory | `transition_matrix()`, `pseudotime_from_transition()`, `build_state_features()` |
| Single-cell | `runIOT()` — matrix / SingleCellExperiment / Seurat |
| Visualisation | `plot_transition_heatmap()`, `plot_transition_flow()`, `plot_theta()`, `plot_pathway_trend()` |
| Demo data | `simulate_iot_states()`, `demo_iot_states` |
| Bulk cohorts | `pathway_markers`, `score_pathways()`, `collapse_probes()`, `gsm_id()` |

## How it works

bioIOT solves

```text
min_P  <C, P> − eps·H(P) + mu·KL(col(P) ‖ b)    s.t.  P·1 = a
```

with a hard source-side row marginal and a KL-anchored column marginal:

- `mu → ∞` recovers hard-marginal OT (pure column features unidentifiable);
- `mu → 0` recovers a plain row-softmax (no target-composition anchoring);
- finite `mu` interpolates the two — the paper's working point is
  `mu = 0.5, eps = 1.0, lam = 0.05`.

A vignette is bundled with the package
(`browseVignettes("bioIOT")` after installation).

## Testing

```r
library(testthat); library(bioIOT)
test_dir("tests/testthat")   # from the repository root
```

## Citation

If you use bioIOT, please cite:

```bibtex
@misc{dong2026bioiotr,
  author       = {Dong, Han},
  title        = {bioIOT: Inverse Optimal Transport for Single-Cell
                  Trajectory Analysis},
  year         = {2026},
  howpublished = {\url{https://github.com/XTSgreen/bioIOT-R}},
  note         = {R package version 0.2.0}
}
```

## License

[MIT](LICENSE) © 2026 Han Dong (XTSgreen)
