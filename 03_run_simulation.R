## =====================================================================
## 03_run_simulation.R
## Full simulation: every scenario in the table below, R replicates each.
## Each scenario is saved to its own file as soon as it finishes, together
## with its settings.  On a restart, a saved scenario is skipped ONLY if its
## saved settings (including R and the seed) match the current ones;
## otherwise it is run again.  So the script can be stopped and restarted
## (e.g. if the computer sleeps) without losing work or mixing old results.
## Output: results/sim/scenario_XX.rds, results/scenarios.csv
## =====================================================================

source("R/methods.R")
source("R/simulate.R")
dir.create("results/sim", recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("metadat", quietly = TRUE) &&
    "dat.nissen2007" %in% data(package = "metadat")$results[, "Item"]) {
  d <- metadat::dat.nissen2007
} else {
  load("data/dat.nissen2007.rda"); d <- dat.nissen2007
}

## ---- Settings -------------------------------------------------------------
R     <- 5000      # replicates per scenario (try R <- 20 first as a test run)
CORES <- max(1, parallel::detectCores() - 2)

rate_mi  <- sum(d$cont.infarction) / sum(d$cont.total)   # about 0.0059
rate_cvd <- sum(d$cont.death)      / sum(d$cont.total)   # about 0.0018

## ---- Scenarios ------------------------------------------------------------
## Columns:
##   outcome : event rate to mimic (heart attack or cardiovascular death)
##   OR      : true odds ratio
##   p0_sd   : between-trial spread of control risk (logit scale); 0 = same
##             risk in every trial, 1 = trials differ, as in real data
##   sparsity: multiplier on the event rate (1 = real rate, 0.25 = 4x rarer)
##   alloc   : "real" = real group sizes; "balanced" = each trial split 1:1
##
## Part A (core): both outcomes x three ORs, realistic settings.
core <- expand.grid(outcome = c("heart_attack", "cv_death"),
                    OR_level = c("null", "published", "large"),
                    p0_sd = 1, sparsity = 1, alloc = "real",
                    stringsAsFactors = FALSE)
core$OR <- with(core, ifelse(OR_level == "null", 1,
                       ifelse(OR_level == "large", 2,
                       ifelse(outcome == "heart_attack", 1.4, 1.6))))

## Part B (one factor at a time, heart-attack rate, OR = 1 and 2):
##   same baseline risk in every trial; sparser data; balanced allocation.
vary <- rbind(
  data.frame(p0_sd = 0, sparsity = 1,    alloc = "real"),
  data.frame(p0_sd = 1, sparsity = 0.5,  alloc = "real"),
  data.frame(p0_sd = 1, sparsity = 0.25, alloc = "real"),
  data.frame(p0_sd = 1, sparsity = 1,    alloc = "balanced"))
vary <- merge(vary, data.frame(OR = c(1, 2)))
vary$outcome  <- "heart_attack"
vary$OR_level <- ifelse(vary$OR == 1, "null", "large")

## Part C (added after the first run): extreme sparsity, to find where the
##   methods - in particular unconditional vs conditional Firth - start to
##   differ.  Heart-attack rate 10x and 20x rarer, OR = 1 and 2.
##   Added at the END so scenarios 01-14 keep their IDs and seeds, and their
##   saved results are reused.
extreme <- merge(data.frame(p0_sd = 1, sparsity = c(0.1, 0.05), alloc = "real"),
                 data.frame(OR = c(1, 2)))
extreme$outcome  <- "heart_attack"
extreme$OR_level <- ifelse(extreme$OR == 1, "null", "large")

cols <- c("outcome", "OR_level", "OR", "p0_sd", "sparsity", "alloc")
scenarios <- rbind(core[, cols], vary[, cols], extreme[, cols])
scenarios$id   <- sprintf("%02d", seq_len(nrow(scenarios)))
scenarios$seed <- 1000 + seq_len(nrow(scenarios))
scenarios$R    <- R
write.csv(scenarios, "results/scenarios.csv", row.names = FALSE)
cat("Scenarios to run:\n"); print(scenarios, row.names = FALSE)

## ---- Run ----------------------------------------------------------------------
for (i in seq_len(nrow(scenarios))) {
  sc   <- scenarios[i, ]
  file <- sprintf("results/sim/scenario_%s.rds", sc$id)
  settings <- sc[, c("outcome", "OR", "p0_sd", "sparsity", "alloc", "seed", "R")]
  if (file.exists(file)) {
    old <- attr(readRDS(file), "settings")
    if (!is.null(old) && isTRUE(all.equal(old, settings, check.attributes = FALSE))) {
      cat("\nScenario", sc$id, "already done with the same settings - skipping\n"); next
    }
    cat("\nScenario", sc$id, "was saved with different settings - running again")
  }

  rate <- if (sc$outcome == "heart_attack") rate_mi else rate_cvd
  if (sc$alloc == "balanced") {
    tot <- d$treat.total + d$cont.total
    n1  <- ceiling(tot / 2); n0 <- floor(tot / 2)
  } else {
    n1 <- d$treat.total; n0 <- d$cont.total
  }

  cat(sprintf("\nScenario %s: %s, OR = %.1f, p0_sd = %g, sparsity = %g, alloc = %s ... ",
              sc$id, sc$outcome, sc$OR, sc$p0_sd, sc$sparsity, sc$alloc))
  t0 <- Sys.time()
  sim <- run_sim(R = sc$R, OR = sc$OR, p0 = rate * sc$sparsity,
                 n1 = n1, n0 = n0, seed = sc$seed, cores = CORES,
                 p0_sd = sc$p0_sd)
  sim$scenario <- sc$id
  attr(sim, "settings") <- settings
  saveRDS(sim, file)
  cat(sprintf("done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("\n\nAll scenarios finished. Results are in results/sim/\n")
