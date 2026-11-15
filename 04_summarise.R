## =====================================================================
## 04_summarise.R
## Summarise every scenario saved by 03_run_simulation.R.
## Output: results/summary_all.csv     (one row per scenario x method)
##         results/firth_comparison.csv (unconditional vs conditional Firth)
## =====================================================================

source("R/simulate.R")   # for summarise_sim()

scenarios <- read.csv("results/scenarios.csv", colClasses = c(id = "character"))
files     <- sprintf("results/sim/scenario_%s.rds", scenarios$id)

## A saved scenario is used only if it exists AND was produced with exactly
## the settings currently listed in results/scenarios.csv (outcome, OR,
## p0_sd, sparsity, alloc, seed, R).  Anything else - missing, from an older
## run, or saved before settings were recorded - is left out and listed.
setting_cols <- c("outcome", "OR", "p0_sd", "sparsity", "alloc", "seed", "R")
status <- vapply(seq_len(nrow(scenarios)), function(i) {
  if (!file.exists(files[i])) return("missing")
  old <- attr(readRDS(files[i]), "settings")
  if (is.null(old)) return("no saved settings")
  if (!isTRUE(all.equal(old[, setting_cols], scenarios[i, setting_cols],
                        check.attributes = FALSE))) return("settings differ")
  "ok"
}, character(1))
done <- status == "ok"
if (!all(done)) {
  cat("WARNING: these scenarios are NOT included in the summary",
      "(rerun 03_run_simulation.R to complete them):\n")
  print(data.frame(id = scenarios$id[!done], status = status[!done]), row.names = FALSE)
}
if (!any(done)) stop("No valid scenario results to summarise.")

## ---- 1. Main summary: bias, coverage, rejection rate, failures ----------
## All rates are computed among successful fits (n_ok of R replicates);
## n_failed and n_warning are reported alongside so the denominator is clear.
summary_all <- do.call(rbind, lapply(which(done), function(i) {
  sim <- readRDS(files[i])
  cbind(scenarios[i, c("id", "outcome", "OR", "p0_sd", "sparsity", "alloc")],
        summarise_sim(sim),
        avg_zero_any    = mean(sim$zero_any),
        avg_zero_double = mean(sim$zero_double), row.names = NULL)
}))
summary_all$reject_label <- ifelse(summary_all$OR == 1, "false_positive", "power")
write.csv(summary_all, "results/summary_all.csv", row.names = FALSE)

show <- function(tab, cols) print(tab[, cols], digits = 3, row.names = FALSE)

cat("\n===== FALSE-POSITIVE RATE (true OR = 1; target 0.05) =====\n")
show(subset(summary_all, OR == 1),
     c("id", "outcome", "p0_sd", "sparsity", "alloc", "method",
       "reject", "reject_mcse", "coverage", "n_ok", "n_failed"))

cat("\n===== POWER AND COVERAGE (true OR > 1; coverage target 0.95) =====\n")
show(subset(summary_all, OR > 1),
     c("id", "outcome", "OR", "p0_sd", "sparsity", "alloc", "method",
       "reject", "reject_mcse", "coverage", "cover_mcse", "n_ok", "n_failed"))

cat("\n===== BIAS IN THE LOG ODDS RATIO =====\n")
show(summary_all,
     c("id", "outcome", "OR", "sparsity", "alloc", "method", "bias", "bias_mcse"))

## ---- 2. Power at the real rosiglitazone trial sizes -----------------------
## Scenarios 03 (heart attack, OR 1.4) and 04 (CV death, OR 1.6) use the real
## trial sizes, real event rates and the published effect sizes.
cat("\n===== POWER AT THE REAL TRIAL SIZES (published effect sizes) =====\n")
show(subset(summary_all, alloc == "real" & sparsity == 1 & p0_sd == 1 &
              OR %in% c(1.4, 1.6)),
     c("outcome", "OR", "method", "reject", "reject_mcse"))

## ---- 3. Unconditional vs conditional Firth, paired --------------------------
## Both are fitted to the same simulated data sets, so their difference is
## estimated much more precisely than comparing the two summaries separately.
firth_cmp <- do.call(rbind, lapply(which(done), function(i) {
  sim <- readRDS(files[i])
  u <- sim[sim$method == "UFirth", ]; cf <- sim[sim$method == "CFirth", ]
  u <- u[order(u$rep), ]; cf <- cf[order(cf$rep), ]
  stopifnot(identical(u$rep, cf$rep))
  ok <- !u$error & !cf$error                         # both fits valid
  dd <- u$logOR[ok] - cf$logOR[ok]                  # paired difference
  w  <- (u$hi - u$lo)[ok] - (cf$hi - cf$lo)[ok]     # difference in CI width
  cbind(scenarios[i, c("id", "outcome", "OR", "p0_sd", "sparsity", "alloc")],
        n_pairs = sum(ok),
        mean_diff_logOR  = mean(dd), mcse_diff_logOR  = sd(dd) / sqrt(sum(ok)),
        mean_diff_CIwidth = mean(w), mcse_diff_CIwidth = sd(w) / sqrt(sum(ok)),
        row.names = NULL)
}))
cat("\n===== UNCONDITIONAL minus CONDITIONAL FIRTH (paired, same data sets) =====\n")
cat("Positive mean_diff_logOR = unconditional Firth gives larger estimates\n")
print(firth_cmp, digits = 3, row.names = FALSE)
write.csv(firth_cmp, "results/firth_comparison.csv", row.names = FALSE)

cat("\nSaved to results/summary_all.csv and results/firth_comparison.csv\n")
