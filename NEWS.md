# bioIOT 0.2.0

* New R core solver: `soft_sinkhorn()` (Anderson-accelerated semi-relaxed
  soft Sinkhorn), `row_conditional()`, `make_cost()`, `row_ce_loss()`,
  `zscore_phi()`.
* New fitting: `fit_iot()` with exact implicit gradients, two-stage debias
  refit and multi-restart; `print`/`summary` methods for `bioIOT_fit`.
* New trajectory layer: `transition_matrix()`,
  `pseudotime_from_transition()`, `build_state_features()`.
* New single-cell interfaces: `runIOT()` generic with matrix,
  SingleCellExperiment and Seurat methods.
* New visualisation (ggplot2): `plot_transition_heatmap()`,
  `plot_transition_flow()`, `plot_theta()`, `plot_pathway_trend()`.
* New reproducible demo: `simulate_iot_states()` plus bundled dataset
  `demo_iot_states`.
* New vignette and expanded test suite (implicit-gradient vs finite
  differences, recovery, interfaces, plots).

# bioIOT 0.1.0

* Initial release: `pathway_markers`, `score_pathways`, `collapse_probes`,
  `gsm_id`, `has_cel_file`, `find_platform_file`.
