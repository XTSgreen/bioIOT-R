# Bundled demo dataset: full pipeline + visualisation
if (nzchar(Sys.getenv("BIOIOT_LIB"))) {
  base::.libPaths(c(Sys.getenv("BIOIOT_LIB"), base::.libPaths()))
}
suppressPackageStartupMessages({library(bioIOT); library(ggplot2)})
outdir <- "demo_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

data(demo_iot_states)
d <- demo_iot_states
cat(sprintf("[demo] K=%d states, F=%d features, %d cells\n",
            d$K, ncol(d$phi), nrow(d$cell_embedding)))

## fit
fit <- fit_iot(d$phi, d$a, d$b, d$T_true, n_restart = 4, epochs = 300, seed = 1)
Q <- transition_matrix(fit)
pt <- pseudotime_from_transition(Q, root = "S1")
pt_true <- pseudotime_from_transition(d$T_true, root = "S1")

cat("\n=== theta recovery ===\n")
print(data.frame(feature = c("f1", "f2", "sim"),
                 fitted = round(fit$theta, 3),
                 truth = d$theta_true))
cat(sprintf("Pearson corr(theta, theta_true): %.4f\n",
            cor(fit$theta, d$theta_true)))
cat(sprintf("mean per-row L1(Q, T_true): %.4f\n", mean(abs(Q - d$T_true))))
cat(sprintf("pseudotime corr (Spearman, fitted vs true): %.4f\n",
            cor(unlist(pt), unlist(pt_true), method = "spearman")))

## fig1: transition heatmap (fitted)
ggsave(file.path(outdir, "demo_fig1_heatmap.png"),
       plot_transition_heatmap(Q) +
         ggtitle("demo_iot_states: fitted transitions"),
       width = 5.4, height = 4.6, dpi = 150, bg = "white")

## fig2: truth heatmap for comparison
ggsave(file.path(outdir, "demo_fig2_true_heatmap.png"),
       plot_transition_heatmap(d$T_true) +
         ggtitle("demo_iot_states: true transitions"),
       width = 5.4, height = 4.6, dpi = 150, bg = "white")

## fig3: flow arrows on the 2-D state embedding
ggsave(file.path(outdir, "demo_fig3_flow.png"),
       plot_transition_flow(Q, d$embedding, threshold = 0.08) +
         ggtitle("demo_iot_states: IOT flow"),
       width = 5.4, height = 4.6, dpi = 150, bg = "white")

## fig4: feature weights vs truth
ggsave(file.path(outdir, "demo_fig4_theta.png"),
       plot_theta(fit, labels = c("f1", "f2", "sim")) +
         ggtitle("demo_iot_states: fitted weights (truth: 0.9 / -0.7 / 1.1)"),
       width = 5.0, height = 3.8, dpi = 150, bg = "white")

## fig5: pseudotime, fitted vs true
df <- data.frame(state = names(pt_true),
                 truth = unlist(pt_true),
                 fitted = unlist(pt[names(pt_true)]))
p5 <- ggplot(df, aes(truth, fitted, label = state)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(color = "grey55") +
  geom_point(color = "#B2182B", size = 3.4) +
  geom_text(vjust = -1.2, size = 3.4) +
  labs(x = "pseudotime under TRUE T", y = "pseudotime under fitted Q",
       title = "demo_iot_states: pseudotime recovery (root = S1)") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "demo_fig5_pseudotime.png"), p5,
       width = 5.2, height = 4.4, dpi = 150, bg = "white")

cat("\n[done] outputs in", normalizePath(outdir), "\n")
