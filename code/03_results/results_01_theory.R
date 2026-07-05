# =============================================================================
# results_01_theory.R
# Author : JoonHo Lee (jlee296@ua.edu)
# -----------------------------------------------------------------------------
# Figure 1: the closed-form Design Effect Ratio curves (Theorems 1-2). These are
# analytic -- no data, no model fitting. The random-effect DER is an inverted-U
# in the shrinkage factor B; the fixed-effect DER spans the band
# [DEFF(1 - B), DEFF], with within-group covariates at the upper edge and
# between-group covariates at the lower edge.
#
# Writes : output/figures/fig1_theory.pdf
# Usage  : Rscript code/03_results/results_01_theory.R
# =============================================================================
suppressMessages(library(ggplot2))

ROOT <- tryCatch(here::here(), error = function(e) getwd())
out  <- file.path(ROOT, "output", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
J <- 50
kappa <- function(B, J) 1 - 1 / (J * (1 - B) + B)
Bg <- seq(0.001, 0.999, length.out = 400)
deffs <- c(1.5, 2, 3, 5)

# Left panel: random-effect DER = B * DEFF * kappa_J(B), + large-J limit
dl <- do.call(rbind, lapply(deffs, function(D) {
  rbind(
    data.frame(B = Bg, der = Bg * D * kappa(Bg, J), DEFF = D, type = "exact (J = 50)"),
    data.frame(B = Bg, der = Bg * D, DEFF = D, type = "large-J limit")
  )
}))
dl$DEFF <- factor(dl$DEFF, labels = paste0("DEFF = ", deffs))

pL <- ggplot(dl, aes(B, der, color = DEFF, linetype = type)) +
  geom_hline(yintercept = c(1.2, 1.5), color = "gray55", linewidth = 0.35) +
  annotate("text", x = 0.03, y = 1.2, label = "c[0] == 1.2", parse = TRUE,
           hjust = 0, vjust = -0.5, size = 3.1, color = "gray30") +
  annotate("text", x = 0.03, y = 1.5, label = "c[0] == 1.5", parse = TRUE,
           hjust = 0, vjust = -0.5, size = 3.1, color = "gray30") +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = okabe) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  labs(x = expression(paste("shrinkage factor ", B[j])),
       y = expression(DER[j]),
       title = "Random effects: DER = B %.% DEFF %*% kappa(J)" ) +
  labs(title = expression(paste("Random effects:  ", DER[j] == B[j] %.% DEFF[j] %.% kappa[j](J)))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.box = "vertical",
        plot.title = element_text(size = 11))

# Right panel: fixed-effect band [DEFF(1-B), DEFF]
dr <- do.call(rbind, lapply(c(2, 5), function(D) {
  data.frame(B = Bg, lower = D * (1 - Bg), upper = D, DEFF = D)
}))
dr$DEFF <- factor(dr$DEFF, labels = paste0("DEFF = ", c(2, 5)))

pR <- ggplot(dr, aes(B)) +
  geom_hline(yintercept = c(1.2, 1.5), color = "gray55", linewidth = 0.35) +
  annotate("text", x = 0.03, y = 1.2, label = "c[0] == 1.2", parse = TRUE,
           hjust = 0, vjust = -0.5, size = 3.1, color = "gray30") +
  annotate("text", x = 0.03, y = 1.5, label = "c[0] == 1.5", parse = TRUE,
           hjust = 0, vjust = -0.5, size = 3.1, color = "gray30") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = DEFF), alpha = 0.25) +
  geom_line(aes(y = upper, color = DEFF), linewidth = 0.7) +
  geom_line(aes(y = lower, color = DEFF), linewidth = 0.7, linetype = "22") +
  scale_color_manual(values = okabe[c(2, 4)]) +
  scale_fill_manual(values = okabe[c(2, 4)]) +
  labs(x = expression(paste("shrinkage factor ", B)),
       y = expression(DER[beta[k]]),
       title = expression(paste("Fixed effects:  ", DER[beta[k]] == DEFF %.% (1 - R[k]), ",  ", R[k] %in% "[0, B]"))) +
  annotate("text", x = 0.55, y = 4.6, label = "within-group (solid)",
           size = 3.1, color = "gray25") +
  annotate("text", x = 0.72, y = 0.6, label = "between-group (dashed)",
           size = 3.1, color = "gray25") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 11))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p <- pL + pR
  ggsave(file.path(out, "fig1_theory.pdf"), p, width = 10, height = 4.6)
  ggsave(file.path(out, "fig1_theory.png"), p, width = 10, height = 4.6, dpi = 150)
} else {
  pdf(file.path(out, "fig1_theory.pdf"), width = 10, height = 4.6)
  gridExtra::grid.arrange(pL, pR, ncol = 2)
  dev.off()
}
cat("written:", file.path(out, "fig1_theory.pdf"), "\n")
