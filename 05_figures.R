## =====================================================================
## 05_figures.R
## Figures and a LaTeX table for the paper (base R graphics, no extra packages).
## Run after 01_real_data.R and 04_summarise.R.
## Output (results/figures/):
##   fig1_real_data.pdf       odds ratio and 95% CI from every method, real data
##   fig2_false_positive.pdf  false-positive rate by scenario (true OR = 1)
##   fig3_coverage.pdf        CI coverage in every scenario
##   fig4_sparsity_bias.pdf   bias as events become rarer
##   table_real_data.tex      LaTeX table of the real-data results
## Every simulation estimate is drawn with a 95% Monte Carlo interval
## (estimate +/- 1.96 x Monte Carlo standard error).
## =====================================================================

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

real <- read.csv("results/real_data_results.csv")
summ <- read.csv("results/summary_all.csv", colClasses = c(id = "character"))

## ---- Shared style ---------------------------------------------------------
## Colour-blind-safe palette (validated), plus a distinct marker shape per
## method so the figures also read in grayscale.
meth_order <- c("MH", "MH05", "Peto", "UFirth", "CFirth")
meth_label <- c(MH = "Mantel-Haenszel", MH05 = "MH + 0.5 correction",
                Peto = "Peto", UFirth = "Unconditional Firth",
                CFirth = "Conditional Firth")
meth_col   <- c(MH = "#2a78d6", MH05 = "#eb6834", Peto = "#1baf7a",
                UFirth = "#eda100", CFirth = "#e87ba4")
meth_pch   <- c(MH = 16, MH05 = 17, Peto = 15, UFirth = 18, CFirth = 1)
ink   <- "#2b2b2b"; muted <- "#8a8a85"; grid_col <- "#e6e6e3"

setup <- function() par(mar = c(4.2, 9.5, 2.2, 1), mgp = c(2.4, 0.6, 0),
                        tcl = -0.3, col.axis = ink, col.lab = ink,
                        family = "sans", cex.axis = 0.85, las = 1)

## Legend in its own strip below the plot(s), so it never covers data
legend_strip <- function(ncol = 5) {
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = meth_label[meth_order], col = meth_col[meth_order],
         pch = meth_pch[meth_order], lwd = 1.5, bty = "n", cex = 0.8,
         text.col = ink, ncol = ncol)
}

scen_label <- function(s) {
  out <- ifelse(s$outcome == "heart_attack", "MI", "CVD")
  out <- paste0(out, " OR ", s$OR)
  extra <- character(nrow(s))
  extra[s$p0_sd == 0]         <- "same risk"
  sp <- s$sparsity < 1
  extra[sp] <- paste0(round(1 / s$sparsity[sp]), "x rarer")
  extra[s$alloc == "balanced"] <- "1:1 groups"
  ifelse(extra == "", out, paste0(out, ", ", extra))
}

## Dot plot with Monte Carlo error bars: one row per scenario, one dot per method
dotplot_mc <- function(tab, value, mcse, ref, xlab, main) {
  tab <- tab[order(tab$id), ]
  ids <- unique(tab$id)
  labs <- scen_label(tab[match(ids, tab$id), ])
  y0 <- rev(seq_along(ids))                         # first scenario at top
  off <- seq(-0.28, 0.28, length.out = length(meth_order))
  tab <- tab[is.finite(tab[[value]]), ]
  lo <- tab[[value]] - 1.96 * tab[[mcse]]; hi <- tab[[value]] + 1.96 * tab[[mcse]]
  setup()
  plot(NA, xlim = range(c(lo, hi, ref), na.rm = TRUE), ylim = c(0.4, length(ids) + 0.6),
       yaxt = "n", xlab = xlab, ylab = "", main = main, cex.main = 0.95,
       font.main = 1, bty = "n")
  abline(h = y0, col = grid_col, lwd = 0.8)
  abline(v = ref, lty = 2, col = muted)
  axis(2, at = y0, labels = labs, tick = FALSE, cex.axis = 0.8)
  for (j in seq_along(meth_order)) {
    m <- meth_order[j]; s <- tab[tab$method == m, ]
    y <- y0[match(s$id, ids)] + off[j]
    segments(s[[value]] - 1.96 * s[[mcse]], y, s[[value]] + 1.96 * s[[mcse]], y,
             col = meth_col[m], lwd = 1.5)
    points(s[[value]], y, pch = meth_pch[m], col = meth_col[m], cex = 1.1)
  }
}

## ---- Figure 1: real data ------------------------------------------------------
pdf("results/figures/fig1_real_data.pdf", width = 7, height = 3.6)
par(mfrow = c(1, 2))
for (o in c("heart_attack", "cv_death")) {
  r <- real[real$outcome == o, ]; r <- r[match(rev(meth_order), r$method), ]
  setup(); par(mar = c(4.2, 9.5, 2.2, 0.8))
  xr <- range(c(1, r$CI_lower, r$CI_upper), finite = TRUE) * c(0.95, 1.05)
  plot(NA, xlim = xr, ylim = c(0.5, nrow(r) + 0.5), log = "x",
       yaxt = "n", ylab = "", xlab = "Odds ratio (log scale)", bty = "n",
       main = ifelse(o == "heart_attack", "Heart attack", "Cardiovascular death"),
       cex.main = 0.95, font.main = 1)
  abline(v = 1, lty = 2, col = muted)
  abline(h = seq_len(nrow(r)), col = grid_col, lwd = 0.8)
  axis(2, at = seq_len(nrow(r)), labels = meth_label[r$method], tick = FALSE, cex.axis = 0.8)
  segments(r$CI_lower, seq_len(nrow(r)), r$CI_upper, seq_len(nrow(r)),
           col = meth_col[r$method], lwd = 2)
  points(r$OR, seq_len(nrow(r)), pch = meth_pch[r$method], col = meth_col[r$method], cex = 1.2)
}
dev.off()

## ---- Figure 2: false-positive rate --------------------------------------------
pdf("results/figures/fig2_false_positive.pdf", width = 7, height = 5.6)
layout(matrix(1:2, 2), heights = c(1, 0.12))
dotplot_mc(subset(summ, OR == 1), "reject", "reject_mcse", ref = 0.05,
           xlab = "False-positive rate (target 0.05)",
           main = "Rejection rate when the drug has no effect (true OR = 1)")
legend_strip(ncol = 3)
dev.off()

## ---- Figure 3: coverage ---------------------------------------------------------
pdf("results/figures/fig3_coverage.pdf", width = 7, height = 7.5)
layout(matrix(1:2, 2), heights = c(1, 0.09))
dotplot_mc(summ, "coverage", "cover_mcse", ref = 0.95,
           xlab = "Coverage of the 95% confidence interval (target 0.95)",
           main = "Confidence interval coverage, all scenarios")
legend_strip(ncol = 3)
dev.off()

## ---- Figure 4: bias as events become rarer -----------------------------------
sp <- subset(summ, outcome == "heart_attack" & alloc == "real" & p0_sd == 1)
pdf("results/figures/fig4_sparsity_bias.pdf", width = 7, height = 4.1)
layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(1, 0.16))
for (or in c(1, 2)) {
  s <- sp[sp$OR == or, ]
  setup(); par(mar = c(4.2, 4.2, 2.2, 0.8))
  xs <- sort(unique(s$sparsity))
  yl <- range(c(s$bias - 1.96 * s$bias_mcse, s$bias + 1.96 * s$bias_mcse, 0),
              finite = TRUE)
  plot(NA, xlim = range(xs), ylim = yl, log = "x", xaxt = "n", bty = "n",
       xlab = "Event rate (relative to real data)", ylab = "Bias in log OR",
       main = sprintf("True OR = %g", or), cex.main = 0.95, font.main = 1)
  axis(1, at = xs, labels = paste0(xs, "x"))
  abline(h = 0, lty = 2, col = muted)
  for (m in meth_order) {
    t <- s[s$method == m & is.finite(s$bias), ]; t <- t[order(t$sparsity), ]
    if (nrow(t) == 0) next
    lines(t$sparsity, t$bias, col = meth_col[m], lwd = 2)
    arrows(t$sparsity, t$bias - 1.96 * t$bias_mcse, t$sparsity, t$bias + 1.96 * t$bias_mcse,
           angle = 90, code = 3, length = 0.03, col = meth_col[m])
    points(t$sparsity, t$bias, pch = meth_pch[m], col = meth_col[m], cex = 1.1)
  }
}
legend_strip(ncol = 3)
dev.off()

## ---- LaTeX table: real-data results ---------------------------------------------
fmt <- function(r) sprintf("%.2f (%.2f, %.2f) & %.3f", r$OR, r$CI_lower, r$CI_upper, r$p)
rows <- sapply(meth_order, function(m) {
  a <- real[real$outcome == "heart_attack" & real$method == m, ]
  b <- real[real$outcome == "cv_death"     & real$method == m, ]
  sprintf("%s & %s & %s \\\\", meth_label[m], fmt(a), fmt(b))
})
tex <- c("\\begin{table}[t]",
         "\\caption{Odds ratios (95\\% CI) and $p$-values from each method applied to the rosiglitazone trials.}",
         "\\label{tab:realdata}",
         "\\centering\\small",
         "\\begin{tabular}{lcccc}",
         "\\toprule",
         " & \\multicolumn{2}{c}{Heart attack} & \\multicolumn{2}{c}{Cardiovascular death} \\\\",
         "Method & OR (95\\% CI) & $p$ & OR (95\\% CI) & $p$ \\\\",
         "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, "results/figures/table_real_data.tex")

cat("Saved figures and table to results/figures/\n")
