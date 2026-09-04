# bioIOT showcase: end-to-end run on a real public dataset
#
# Data    : HSMMSingleCell (Bioconductor) -- human skeletal muscle myoblast
#           differentiation time course, 271 cells x ~18k genes (FPKM),
#           timepoints 0 / 24 / 48 / 72 h (Shin et al. 2015; Artistic-2.0,
#           see LICENSE_AUDIT.md)
# Pipeline: log2-FPKM -> HVG -> PCA -> state clustering -> nearest-centroid
#           T_obs per consecutive timepoint pair -> multi-scenario fit_iot
#           -> transition matrices -> random-walk pseudotime
# Bench   : slingshot pseudotime (field standard) vs bioIOT pseudotime,
#           both against the known Hours labels (Spearman)
# Run     : BIOIOT_LIB=<lib> Rscript run_hsmm.R

if (nzchar(Sys.getenv("BIOIOT_LIB"))) {
  base::.libPaths(c(Sys.getenv("BIOIOT_LIB"), base::.libPaths()))
}
suppressPackageStartupMessages({
  library(bioIOT)
  library(HSMMSingleCell)
  library(ggplot2)
})
outdir <- "hsmm_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------- data ----
data("HSMM_expr_matrix")
data("HSMM_sample_sheet")
expr <- as.matrix(HSMM_expr_matrix)
sheet <- as.data.frame(HSMM_sample_sheet)
cells <- intersect(colnames(expr), rownames(sheet))
expr <- expr[, cells]
sheet <- sheet[cells, ]
hours <- as.numeric(as.character(sheet$Hours))  # factor "0"/"24"/"48"/"72"
cat(sprintf("[data] %d genes x %d cells; timepoints: %s\n",
            nrow(expr), ncol(expr), paste(sort(unique(hours)), collapse = "/")))

## ------------------------------------------------------------ features ----
log_expr <- log2(expr + 1)
keep <- rowMeans(log_expr) > 0.05 & apply(log_expr, 1, sd) > 0.4
log_expr <- log_expr[keep, ]
hvg <- head(order(apply(log_expr, 1, var), decreasing = TRUE), 1000)
pca <- prcomp(t(log_expr[hvg, ]), center = TRUE, scale. = FALSE)
n_pc <- min(8, sum(pca$sdev > 1e-8))
emb_all <- pca$x[, seq_len(n_pc)]
colnames(emb_all) <- paste0("PC", seq_len(n_pc))
cat(sprintf("[pca] %d HVGs, PC1-2 var share %.1f%%\n",
            length(hvg), 100 * summary(pca)$importance[3, 2]))

## -------------------------------------------------------------- states ----
set.seed(1)
km <- kmeans(emb_all, centers = 6, nstart = 25)
cl <- km$cluster
ord <- order(tapply(hours, cl, mean))            # relabel so S1 = earliest
lab <- setNames(sprintf("S%d", seq_along(ord)), ord)
state <- factor(lab[as.character(cl)], levels = sprintf("S%d", seq_along(ord)))
cat("[states] per-state n:\n")
print(table(state))
# state centroids over ALL cells (state definition is timepoint-independent;
# empty source timepoints simply give a_i = 0, which bioIOT masks legally)
u_all <- t(vapply(levels(state), function(s)
  colMeans(emb_all[state == s, , drop = FALSE]), numeric(ncol(emb_all))))
rownames(u_all) <- levels(state)

## ------------------------------------------- observed transitions (T_hat) --
# lineage-free heuristic (documented): each target-timepoint cell inherits
# the state of its nearest source-timepoint CELL in PCA space; the resulting
# (source state x target state) counts, with 0.5 additive smoothing against
# empty rows, are row-normalized into a row-conditional T_hat per pair.
pairs <- list(c(0, 24), c(24, 48), c(48, 72))
scen <- lapply(pairs, function(pr) {
  src <- which(hours == pr[1])
  tgt <- which(hours == pr[2])
  src_pos <- src[state[src] %in% levels(state)]
  d <- as.matrix(dist(rbind(emb_all[src_pos, , drop = FALSE],
                            emb_all[tgt, , drop = FALSE])))
  nn <- apply(d[(length(src_pos) + 1):nrow(d), seq_len(length(src_pos)),
                drop = FALSE], 1, which.min)
  from_state <- state[src_pos][nn]         # inherited source state per target cell
  to_state <- state[tgt]
  counts <- table(factor(from_state, levels = levels(state)),
                  factor(to_state, levels = levels(state))) + 0.5
  counts <- as.matrix(counts)
  list(u = u_all,
       a = as.numeric(table(factor(state[src], levels = levels(state)))) / length(src),
       b = as.numeric(table(factor(state[tgt], levels = levels(state)))) / length(tgt),
       T = counts / rowSums(counts))
})

## -------------------------------------------------------------- fitting ----
phi_list <- lapply(scen, function(s) zscore_phi(build_state_features(s$u))$phi_z)
fit <- fit_iot(phi_list,
               lapply(scen, `[[`, "a"),
               lapply(scen, `[[`, "b"),
               lapply(scen, `[[`, "T"),
               mu = 0.5, n_restart = 4, epochs = 300, seed = 1,
               verbose = TRUE)
feat_names <- c(colnames(emb_all), "sim")
cat("\n=== fitted theta (standardized feature space) ===\n")
print(setNames(round(fit$theta, 3), feat_names))
cat("support:", sum(fit$support), "/", fit$F, "features\n")
cat("restart losses:", round(fit$restart_losses, 3), "\n")
write.csv(data.frame(feature = feat_names, theta = fit$theta,
                     support = fit$support),
          file.path(outdir, "hsmm_theta.csv"), row.names = FALSE)

## -------------------------------------------- transitions + pseudotime ----
Qs <- lapply(seq_along(scen), function(i) transition_matrix(fit, which = i))
for (i in seq_along(Qs)) {
  write.csv(Qs[[i]],
            file.path(outdir, sprintf("hsmm_Q_pair%d_%s.csv", i,
                                      paste(pairs[[i]], collapse = "-"))))
}
root <- "S1"
pt_state <- pseudotime_from_transition(Qs[[1]], root = root)
cat("\n=== random-walk pseudotime per state (root =", root, ") ===\n")
print(round(pt_state, 2))
pt_cell <- as.numeric(pt_state[as.character(state)])

## ------------------------------------------------------------ slingshot ----
suppressPackageStartupMessages(library(slingshot))
emb2d <- emb_all[, 1:2]
sds <- slingshot(data.frame(emb2d), clusterLabels = state, start.clus = "S1")
sp <- slingPseudotime(sds)
# HSMM branches into multiple curves; per cell take the mean over its curves
sling_pt <- rowMeans(sp, na.rm = TRUE)
sling_pt[!is.finite(sling_pt)] <- NA
sling_state <- tapply(sling_pt, state, mean, na.rm = TRUE)

comp <- data.frame(
  metric = c("bioIOT pseudotime ~ Hours (cells)",
             "slingshot pseudotime ~ Hours (cells)",
             "bioIOT state pt ~ slingshot state pt"),
  spearman = c(suppressWarnings(cor(pt_cell, hours,
                                    method = "spearman",
                                    use = "complete.obs")),
               suppressWarnings(cor(sling_pt, hours,
                                    method = "spearman",
                                    use = "complete.obs")),
               suppressWarnings(cor(pt_state[levels(state)],
                                    sling_state[levels(state)],
                                    method = "spearman",
                                    use = "complete.obs")))
)
cat("\n=== pseudotime benchmark (Spearman) ===\n")
print(comp, digits = 3)
write.csv(comp, file.path(outdir, "hsmm_pseudotime_benchmark.csv"),
          row.names = FALSE)

## -------------------------------------------------------------- figures ----
ggsave(file.path(outdir, "fig1_transition_heatmap.png"),
       plot_transition_heatmap(Qs[[1]]) +
         ggtitle("HSMM: IOT state transitions, 0h -> 24h"),
       width = 5.8, height = 4.8, dpi = 150)
ggsave(file.path(outdir, "fig2_transition_flow.png"),
       plot_transition_flow(Qs[[1]], scen[[1]]$u[, 1:2], threshold = 0.12) +
         labs(x = "PC1", y = "PC2",
              title = "HSMM: IOT transition flow (0h -> 24h)"),
       width = 5.6, height = 4.8, dpi = 150, bg = "white")
ggsave(file.path(outdir, "fig3_theta.png"),
       plot_theta(fit, labels = feat_names) +
         ggtitle("HSMM: fitted IOT feature weights (3 timepoint pairs)"),
       width = 5.2, height = 4.4, dpi = 150)

df_m <- data.frame(state = levels(state),
                   hours = as.numeric(tapply(hours, state, mean)),
                   bioIOT = as.numeric(pt_state[levels(state)]),
                   slingshot = as.numeric(sling_state[levels(state)]))
p4 <- ggplot(df_m, aes(hours, bioIOT, label = state)) +
  geom_line(color = "grey55") +
  geom_point(color = "#B2182B", size = 3.2) +
  geom_text(vjust = -1.1, size = 3.2) +
  labs(x = "Known time (hours)", y = "bioIOT random-walk pseudotime",
       title = "bioIOT pseudotime tracks the known differentiation clock") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "fig4_pseudotime_vs_hours.png"), p4,
       width = 5.6, height = 4.4, dpi = 150)

myo <- c("MYOD1", "MYOG", "MYH3", "MYF6", "MEF2C", "CKM", "ACTA1", "TNNT3")
# HSMM rownames are Ensembl ids: map to symbols with the shipped annotation
# and collapse to gene level using bioIOT's own collapse_probes()
data("HSMM_gene_annotation")
ann <- as.data.frame(HSMM_gene_annotation)
sym_vec <- ann$gene_short_name[match(rownames(log_expr), rownames(ann))]
ok_g <- !is.na(sym_vec) & nzchar(sym_vec)
log_sym <- collapse_probes(log_expr[ok_g, ], probe_id = rownames(log_expr)[ok_g],
                           gene_id = sym_vec[ok_g])
pw <- score_pathways(log_sym, markers = list(Myogenesis = intersect(myo, rownames(log_sym))))
cat(sprintf("[pathway] myogenic markers found: %d / %d\n",
            sum(myo %in% rownames(log_sym)), length(myo)))
pw_time <- data.frame(hours = hours, score = as.numeric(pw[, "Myogenesis"]))
p5 <- ggplot(pw_time, aes(hours, score)) +
  geom_smooth(method = "loess", se = TRUE, color = "#B2182B",
              fill = "#B2182B22") +
  geom_jitter(width = 1, alpha = 0.45, size = 1) +
  labs(x = "Known time (hours)", y = "Myogenic marker score (z)",
       title = "Myogenic programme over time (pathway util sanity check)") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "fig5_myogenic_trend.png"), p5,
       width = 5.6, height = 4.4, dpi = 150)

cat("\n[done] outputs in", normalizePath(outdir), "\n")
