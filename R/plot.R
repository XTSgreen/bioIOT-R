# bioIOT visualization (ggplot2)

#' Heatmap of a state transition matrix
#'
#' @param Q (K, K) row-stochastic transition matrix (e.g. from
#'   \code{\link{transition_matrix}} or a \code{runIOT} result).
#' @param labels Optional state labels (default: rownames or S1..SK).
#'
#' @return A ggplot object.
#'
#' @examples
#' Q <- matrix(c(0.7, 0.3, 0, 0.2, 0.6, 0.2, 0.1, 0.1, 0.8), 3, 3, byrow = TRUE)
#' plot_transition_heatmap(Q)
#' @export
plot_transition_heatmap <- function(Q, labels = NULL) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  Q <- as.matrix(Q)
  K <- nrow(Q)
  if (is.null(labels)) {
    labels <- if (!is.null(rownames(Q))) rownames(Q) else sprintf("S%d", seq_len(K))
  }
  df <- expand.grid(source = seq_len(K), target = seq_len(K))
  df$weight <- as.vector(Q)  # column-major stacking matches expand.grid order
  df$source <- factor(labels[df$source], levels = labels)
  df$target <- factor(labels[df$target], levels = labels)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$target, y = .data$source,
                                   fill = .data$weight)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$weight)),
                       size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = "#B2182B") +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(x = "Target state", y = "Source state",
                  fill = "Q(i->j)", title = "IOT state transitions") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Flow arrows of state transitions on a 2-D state embedding
#'
#' Draws states as points on a 2-D embedding and transitions as arrows
#' whose width is proportional to the transition mass (CellRank-style).
#'
#' @param Q (K, K) transition matrix.
#' @param embedding (K x >=2) state coordinates (e.g. cluster centroids).
#' @param labels Optional state labels.
#' @param threshold Minimum transition mass to draw an arrow.
#'
#' @return A ggplot object.
#'
#' @examples
#' Q <- matrix(c(0.7, 0.3, 0, 0.2, 0.6, 0.2, 0.1, 0.1, 0.8), 3, 3, byrow = TRUE)
#' emb <- matrix(c(0, 0, 1, 1, 2, 0), 3, 2, byrow = TRUE)
#' plot_transition_flow(Q, emb)
#' @export
plot_transition_flow <- function(Q, embedding, labels = NULL, threshold = 0.05) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  Q <- as.matrix(Q)
  emb <- as.matrix(embedding)
  if (ncol(emb) < 2) stop("`embedding` needs >= 2 columns", call. = FALSE)
  K <- nrow(Q)
  if (is.null(labels)) {
    labels <- if (!is.null(rownames(Q))) rownames(Q) else sprintf("S%d", seq_len(K))
  }
  df <- expand.grid(source = seq_len(K), target = seq_len(K))
  df$w <- as.vector(Q)
  df <- df[df$source != df$target & df$w > threshold, , drop = FALSE]
  pts <- data.frame(x = emb[, 1], y = emb[, 2], label = labels,
                    mass = rowSums(Q))
  p <- ggplot2::ggplot()
  if (nrow(df) > 0) {
    arrow_df <- data.frame(
      x = emb[df$source, 1], y = emb[df$source, 2],
      xend = emb[df$target, 1], yend = emb[df$target, 2], weight = df$w
    )
    p <- p + ggplot2::geom_segment(
      data = arrow_df,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                   yend = .data$yend, linewidth = .data$weight),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.12, "inches"),
                             type = "closed"),
      color = "#2166AC", alpha = 0.85,
      lineend = "round"
    ) +
      ggplot2::scale_linewidth(range = c(0.2, 2.2))
  }
  p + ggplot2::geom_point(data = pts, ggplot2::aes(x = .data$x, y = .data$y,
                                                   size = .data$mass),
                          color = "#B2182B", alpha = 0.9) +
    ggplot2::geom_text(data = pts, ggplot2::aes(x = .data$x, y = .data$y,
                                                label = .data$label),
                       vjust = -1.1, size = 3.4) +
    ggplot2::scale_size(range = c(3, 8), guide = "none") +
    ggplot2::labs(title = "IOT transition flow", linewidth = "Q(i->j)",
                  x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 12)
}

#' Barplot of fitted IOT feature weights
#'
#' @param fit A \code{bioIOT_fit} object (or a named numeric vector of
#'   weights with an optional \code{support} attribute).
#' @param labels Optional feature labels (default: names or f1..fF).
#'
#' @return A ggplot object.
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_iot_states(K = 5, seed = 1)
#' fit <- fit_iot(sim$phi, sim$a, sim$b, sim$T, n_restart = 1, epochs = 100)
#' plot_theta(fit)
#' @export
plot_theta <- function(fit, labels = NULL) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  if (inherits(fit, "bioIOT_fit")) {
    theta <- fit$theta
    support <- fit$support
  } else {
    theta <- as.numeric(fit)
    support <- attr(fit, "support")
    if (is.null(support)) support <- rep(TRUE, length(theta))
  }
  F <- length(theta)
  if (is.null(labels)) {
    labels <- if (!is.null(names(theta))) names(theta) else sprintf("f%d", seq_len(F))
  }
  df <- data.frame(feature = labels, theta = theta,
                   support = as.logical(support), idx = seq_len(F))
  ggplot2::ggplot(df,
                  ggplot2::aes(x = .data$feature, y = .data$theta,
                               fill = .data$support)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_hline(yintercept = 0, color = "grey30") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey70"),
                               guide = "none") +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = expression(theta),
                  title = "IOT feature weights (red = selected support)") +
    ggplot2::theme_minimal(base_size = 12)
}

#' Pathway score trends over time
#'
#' Plots each pathway's score trajectory over a time variable (loess
#' smooth). Works with the output of \code{\link{score_pathways}}.
#'
#' @param pw_scores Sample x pathway score matrix.
#' @param time Numeric time / pseudotime per sample (row).
#' @param se Show the smooth confidence band (default TRUE).
#'
#' @return A ggplot object.
#'
#' @examples
#' set.seed(1)
#' expr <- matrix(rnorm(200 * 12), nrow = 200)
#' rownames(expr) <- paste0("G", 1:200)
#' colnames(expr) <- paste0("S", 1:12)
#' rownames(expr)[1] <- "VIM"; rownames(expr)[2] <- "CDH2"; rownames(expr)[3] <- "MKI67"
#' pw <- score_pathways(expr)
#' plot_pathway_trend(pw, time = sort(runif(12)))
#' @export
plot_pathway_trend <- function(pw_scores, time, se = TRUE) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  pw <- as.matrix(pw_scores)
  time <- as.numeric(time)
  if (length(time) != nrow(pw)) {
    stop("`time` must have one value per sample (row)", call. = FALSE)
  }
  long <- do.call(rbind, lapply(colnames(pw), function(p) {
    data.frame(pathway = p, time = time, score = as.numeric(pw[, p]))
  }))
  long$pathway <- factor(long$pathway, levels = colnames(pw))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$time, y = .data$score,
                                     color = .data$pathway)) +
    ggplot2::geom_smooth(method = "loess", se = se, formula = y ~ x) +
    ggplot2::labs(x = "Time", y = "Pathway score (z)", color = "Pathway",
                  title = "Pathway trajectories") +
    ggplot2::theme_minimal(base_size = 12)
}
