# Synthetic benchmark: can IOT structure denoise noisy observed transitions?
#
# Design (50 replicates, K = 6):
#   truth     T_true from simulate_iot_states (mu = 0.5, eps = 1)
#   noisy obs T_hat  = row-normalized multinomial draws from T_true
#                      (n = 30 cells per source state) -- finite-sampling regime
#   methods   raw      : use T_hat as-is
#             uniform  : 1/K everywhere
#             shrink   : 0.5 * T_hat + 0.5 * uniform
#             quasi-hard : bioIOT fit with mu = 20 (approximates hard-marginal OT)
#             bioIOT   : fit_iot (mu = 0.5, paper working point)
#   metrics   mean per-row L1 distance to T_true; row cross-entropy of T_true
#             under Q_hat; Pearson corr(theta_hat, theta_true) for the IOT fits
# Run: BIOIOT_LIB=<lib> Rscript benchmark_synthetic.R

if (nzchar(Sys.getenv("BIOIOT_LIB"))) {
  base::.libPaths(c(Sys.getenv("BIOIOT_LIB"), base::.libPaths()))
}
library(bioIOT)
library(ggplot2)

set.seed(2026)
outdir <- "benchmark_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# optional CLI override: Rscript benchmark_synthetic.R <n_rep>
args <- commandArgs(trailingOnly = TRUE)
N_REP <- if (length(args) >= 1) as.integer(args[1]) else 50
K <- 6
N_CELLS <- 30

# resume support: already-finished replicates are skipped
raw_path <- file.path(outdir, "benchmark_raw.csv")
done <- if (file.exists(raw_path)) unique(read.csv(raw_path)$rep) else integer(0)
res <- vector("list", N_REP)

for (rep in seq_len(N_REP)) {
  if (rep %in% done) next
  sim <- simulate_iot_states(K = K, seed = rep)
  T_true <- sim$T_true

  # finite-sampling noise: multinomial draws per source state
  T_hat <- t(vapply(seq_len(K), function(i) {
    counts <- rmultinom(1, N_CELLS, prob = T_true[i, ])
    as.numeric(counts) / sum(counts)
  }, numeric(K)))

  # methods
  fit_soft <- fit_iot(sim$phi, sim$a, sim$b, T_hat,
                      mu = 0.5, n_restart = 1, epochs = 200, seed = 1)
  Q_soft <- transition_matrix(fit_soft)
  fit_hard <- fit_iot(sim$phi, sim$a, sim$b, T_hat,
                      mu = 20, n_restart = 1, epochs = 200, seed = 1)
  Q_hard <- transition_matrix(fit_hard)
  Q_raw <- T_hat
  Q_unif <- matrix(1 / K, K, K)
  Q_shrink <- 0.5 * T_hat + 0.5 * Q_unif

  l1 <- function(Q) mean(abs(Q - T_true))
  # Q_hat is already row-conditional (rows sum to 1), so the likelihood of
  # the true transitions under Q_hat is -sum(T_true * log(Q_hat))
  ce <- function(Q) -sum(T_true * log(Q + 1e-12))

  row_i <- data.frame(
    rep = rep,
    method = c("raw", "uniform", "shrink", "quasi-hard (mu=20)", "bioIOT (mu=0.5)"),
    L1 = c(l1(Q_raw), l1(Q_unif), l1(Q_shrink), l1(Q_hard), l1(Q_soft)),
    CE = c(ce(Q_raw), ce(Q_unif), ce(Q_shrink), ce(Q_hard), ce(Q_soft)),
    theta_corr = c(NA, NA, NA,
                   suppressWarnings(cor(fit_hard$theta, sim$theta_true)),
                   suppressWarnings(cor(fit_soft$theta, sim$theta_true))),
    stringsAsFactors = FALSE
  )
  write.table(row_i, raw_path, row.names = FALSE, sep = ",",
              append = file.exists(raw_path),
              col.names = !file.exists(raw_path))
  if (rep %% 10 == 0) cat(sprintf("[bench] %d / %d\n", rep, N_REP))
}

tab <- read.csv(raw_path)
tab <- tab[tab$rep <= N_REP, ]

agg <- do.call(rbind, lapply(unique(tab$method), function(m) {
  s <- tab[tab$method == m, ]
  data.frame(method = m,
             L1_mean = mean(s$L1), L1_sd = sd(s$L1),
             CE_mean = mean(s$CE), CE_sd = sd(s$CE),
             theta_corr_mean = mean(s$theta_corr, na.rm = TRUE),
             theta_corr_sd = sd(s$theta_corr, na.rm = TRUE))
}))
agg <- agg[order(agg$L1_mean), ]
write.csv(agg, file.path(outdir, "benchmark_summary.csv"), row.names = FALSE)
cat(sprintf("\n=== benchmark summary (%d replicates, K=%d, %d cells/state) ===\n",
            length(unique(tab$rep)), K, N_CELLS))
print(agg, digits = 3)

tab$method <- factor(tab$method, levels = agg$method)
p <- ggplot(tab, aes(method, L1, fill = method)) +
  geom_boxplot(outlier.size = 0.8, show.legend = FALSE) +
  coord_flip() +
  labs(x = NULL, y = "mean per-row L1 distance to true T",
       title = sprintf("Transition-matrix recovery under finite-sampling noise (%d reps)",
                       length(unique(tab$rep)))) +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "benchmark_L1.png"), p, width = 7, height = 4, dpi = 150)

p2 <- ggplot(na.omit(tab[, c("method", "theta_corr")]),
             aes(method, theta_corr, fill = method)) +
  geom_boxplot(outlier.size = 0.8, show.legend = FALSE) +
  coord_flip() +
  labs(x = NULL, y = "corr(theta_hat, theta_true) (Pearson)",
       title = "Feature-weight recovery") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "benchmark_theta.png"), p2, width = 7, height = 3, dpi = 150)

cat("\n[done] outputs in", normalizePath(outdir), "\n")
