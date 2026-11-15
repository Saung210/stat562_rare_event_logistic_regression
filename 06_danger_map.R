## =====================================================================
## 06_danger_map.R
## (1) DANGER MAP: across a grid of group imbalance x event rarity, how
##     often does the 0.5 continuity correction miss a real harm that
##     conditional Firth detects?  True OR = 1.4 (the published estimate).
## (2) FRAGILITY INDEX: for the real rosiglitazone data, how many events
##     must change before each method's conclusion flips?
## (3) Rosiglitazone's position is marked on the map.
## Run after 01_real_data.R.  Uses the four fast methods (unconditional
## Firth gives the same estimates as conditional Firth and is much slower).
## Output: results/map/cell_XX.rds, results/map/map_results.csv,
##         results/fragility.csv, results/figures/fig5_danger_map.pdf,
##         results/figures/table_fragility.tex
## =====================================================================

source("R/methods.R")
source("R/simulate.R")
dir.create("results/map", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("metadat", quietly = TRUE) &&
    "dat.nissen2007" %in% data(package = "metadat")$results[, "Item"]) {
  d <- metadat::dat.nissen2007
} else {
  load("data/dat.nissen2007.rda"); d <- dat.nissen2007
}

## ---- Settings -------------------------------------------------------------
R_MAP   <- 2000       # replicates per map cell (try 50 first as a test)
CORES   <- max(1, parallel::detectCores() - 2)
OR_TRUE <- 1.4        # published heart-attack estimate
ratios  <- c(1, 1.5, 2, 3, 4)            # treatment : control allocation
rarity  <- c(1, 0.5, 0.3, 0.2, 0.1)      # event rate relative to real MI rate
MAP_METHODS <- c("MH", "MH05", "Peto", "CFirth")

rate_mi  <- sum(d$cont.infarction) / sum(d$cont.total)
rate_cvd <- sum(d$cont.death)      / sum(d$cont.total)
tot      <- d$treat.total + d$cont.total          # real trial sizes

## =====================================================================
## (1) DANGER MAP
## =====================================================================
grid <- expand.grid(ratio = ratios, rarity = rarity)
grid$cell <- sprintf("%02d", seq_len(nrow(grid)))
grid$seed <- 5000 + seq_len(nrow(grid))

for (i in seq_len(nrow(grid))) {
  g    <- grid[i, ]
  file <- sprintf("results/map/cell_%s.rds", g$cell)
  settings <- data.frame(ratio = g$ratio, rarity = g$rarity, OR = OR_TRUE,
                         seed = g$seed, R = R_MAP)
  if (file.exists(file) &&
      isTRUE(all.equal(attr(readRDS(file), "settings"), settings,
                       check.attributes = FALSE))) {
    cat("Map cell", g$cell, "already done - skipping\n"); next
  }
  ## each trial keeps its real total size, split in the chosen ratio
  n1 <- round(tot * g$ratio / (1 + g$ratio)); n0 <- tot - n1
  cat(sprintf("Map cell %s: ratio %.1f:1, events x%.1f ... ", g$cell, g$ratio, g$rarity))
  t0  <- Sys.time()
  sim <- run_sim(R = R_MAP, OR = OR_TRUE, p0 = rate_mi * g$rarity, n1 = n1, n0 = n0,
                 seed = g$seed, cores = CORES, p0_sd = 1, method_names = MAP_METHODS)
  attr(sim, "settings") <- settings
  saveRDS(sim, file)
  cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## Summarise each cell
map <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  sim <- readRDS(sprintf("results/map/cell_%s.rds", grid$cell[i]))
  w   <- reshape(sim[, c("rep", "method", "p", "logOR", "error")], idvar = "rep",
                 timevar = "method", direction = "wide")
  ok  <- !w$error.MH05 & !w$error.CFirth                 # both methods gave an answer
  sig_f <- w$p.CFirth[ok] < 0.05; sig_c <- w$p.MH05[ok] < 0.05
  hidden <- mean(sig_f & !sig_c)          # Firth detects harm, correction misses it
  pw <- function(m) { e <- w[[paste0("error.", m)]]; mean(w[[paste0("p.", m)]][!e] < 0.05) }
  data.frame(grid[i, c("cell", "ratio", "rarity")], n_pairs = sum(ok),
             hidden_harm = hidden, hidden_mcse = sqrt(hidden * (1 - hidden) / sum(ok)),
             power_MH = pw("MH"), power_MH05 = pw("MH05"),
             power_Peto = pw("Peto"), power_CFirth = pw("CFirth"),
             bias_MH05 = mean(w$logOR.MH05[!w$error.MH05]) - log(OR_TRUE),
             bias_CFirth = mean(w$logOR.CFirth[!w$error.CFirth]) - log(OR_TRUE))
}))
write.csv(map, "results/map/map_results.csv", row.names = FALSE)
cat("\n===== DANGER MAP (true OR = 1.4) =====\n")
cat("hidden_harm = share of meta-analyses where conditional Firth finds the harm\n",
    "             (p < 0.05) but the 0.5 correction does not\n")
print(map[, c("ratio", "rarity", "hidden_harm", "hidden_mcse", "power_CFirth",
              "power_MH05", "bias_MH05")], digits = 3, row.names = FALSE)

## ---- Hidden-harm rate at the REAL rosiglitazone design ----------------------
## The map uses one allocation ratio for every trial; the real trials vary
## (0.5:1 to 4.2:1).  Scenarios 03 and 04 of 03_run_simulation.R use the real
## trial sizes, so they give the hidden-harm rate for the actual design.
real_design <- do.call(rbind, lapply(c("03", "04"), function(id) {
  f <- sprintf("results/sim/scenario_%s.rds", id)
  if (!file.exists(f)) return(NULL)
  s <- readRDS(f)
  w <- reshape(s[, c("rep", "method", "p", "error")], idvar = "rep",
               timevar = "method", direction = "wide")
  ok <- !w$error.MH05 & !w$error.CFirth
  h  <- mean(w$p.CFirth[ok] < 0.05 & w$p.MH05[ok] >= 0.05)
  data.frame(scenario = id,
             outcome = ifelse(id == "03", "heart_attack (OR 1.4)", "cv_death (OR 1.6)"),
             hidden_harm = h, hidden_mcse = sqrt(h * (1 - h) / sum(ok)))
}))
if (!is.null(real_design)) {
  cat("\n===== HIDDEN-HARM RATE AT THE REAL ROSIGLITAZONE DESIGN =====\n")
  print(real_design, digits = 3, row.names = FALSE)
  write.csv(real_design, "results/map/real_design_hidden_harm.csv", row.names = FALSE)
}

## =====================================================================
## (2) FRAGILITY INDEX on the real data
## =====================================================================
## Greedy search (in the spirit of Atal et al., 2019): at each step, change
## ONE event in ONE arm of ONE trial - whichever single change moves the
## p-value furthest towards the 0.05 threshold - until the conclusion flips.
##  * significant result:     remove evidence of harm (one fewer treatment
##                            event, or one more control event)
##  * non-significant result: add evidence of harm (one more treatment
##                            event, or one fewer control event)
## The fragility index is the number of changes needed (capped at 50).
fragility <- function(f, x1, n1, x0, n0, alpha = 0.05, max_steps = 50) {
  pv  <- function(a, c) tryCatch(suppressWarnings(f(a, n1, c, n0)[["p"]]),
                                 error = function(e) NA)
  p0  <- pv(x1, x0); sig0 <- p0 < alpha
  for (step in 1:max_steps) {
    cand <- list()
    for (k in seq_along(x1)) {
      if (sig0) {   # weaken the evidence of harm
        if (x1[k] > 0)     { a <- x1; a[k] <- a[k] - 1; cand[[length(cand) + 1]] <- list(a, x0) }
        if (x0[k] < n0[k]) { c <- x0; c[k] <- c[k] + 1; cand[[length(cand) + 1]] <- list(x1, c) }
      } else {      # strengthen the evidence of harm
        if (x1[k] < n1[k]) { a <- x1; a[k] <- a[k] + 1; cand[[length(cand) + 1]] <- list(a, x0) }
        if (x0[k] > 0)     { c <- x0; c[k] <- c[k] - 1; cand[[length(cand) + 1]] <- list(x1, c) }
      }
    }
    ps <- vapply(cand, function(z) pv(z[[1]], z[[2]]), numeric(1))
    best <- if (sig0) which.max(ps) else which.min(ps)
    x1 <- cand[[best]][[1]]; x0 <- cand[[best]][[2]]
    if ((ps[best] < alpha) != sig0) return(c(fragility = step, p_start = p0, p_end = ps[best]))
  }
  c(fragility = NA, p_start = p0, p_end = NA)
}

frag <- do.call(rbind, lapply(c(heart_attack = "infarction", cv_death = "death"), function(o) {
  x1 <- d[[paste0("treat.", o)]]; x0 <- d[[paste0("cont.", o)]]
  do.call(rbind, lapply(names(methods), function(nm) {
    r <- fragility(methods[[nm]], x1, d$treat.total, x0, d$cont.total)
    data.frame(outcome = ifelse(o == "infarction", "heart_attack", "cv_death"),
               method = nm, p_original = r[["p_start"]],
               significant = r[["p_start"]] < 0.05,
               events_to_flip = r[["fragility"]], p_after = r[["p_end"]])
  }))
}))
write.csv(frag, "results/fragility.csv", row.names = FALSE)
cat("\n===== FRAGILITY INDEX (real rosiglitazone data) =====\n")
cat("events_to_flip = number of single-event changes needed to reverse significance\n")
print(frag, digits = 3, row.names = FALSE)

## =====================================================================
## (3) FIGURE: danger map with rosiglitazone marked
## =====================================================================
ink <- "#2b2b2b"; muted <- "#8a8a85"
ramp <- colorRampPalette(c("#f3f7fd", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                           "#256abf", "#184f95", "#0d366b"))(100)   # one-hue sequential
zmax <- max(map$hidden_harm)
nx <- length(ratios); ny <- length(rarity)

## rosiglitazone's position: overall allocation ratio, event rate relative to MI
rosi <- data.frame(label = c("Rosiglitazone: heart attack", "Rosiglitazone: CV death"),
                   ratio = sum(d$treat.total) / sum(d$cont.total),
                   rarity = c(1, rate_cvd / rate_mi))
pos_x <- function(r) approx(ratios, seq_len(nx), xout = r, rule = 2)$y
pos_y <- function(s) approx(rev(rarity), seq_len(ny), xout = s, rule = 2)$y  # rarer = lower

pdf("results/figures/fig5_danger_map.pdf", width = 7, height = 5.8)
layout(matrix(1:2, 1), widths = c(1, 0.16))
par(mar = c(6.2, 5.2, 3, 1), mgp = c(2.6, 0.6, 0), las = 1, col.axis = ink,
    col.lab = ink, family = "sans", cex.axis = 0.85)
plot(NA, xlim = c(0.5, nx + 0.5), ylim = c(0.5, ny + 0.5), xaxt = "n", yaxt = "n",
     xlab = "Allocation ratio (treatment : control)", ylab = "Event rate (relative to real heart-attack rate)",
     bty = "n", main = "How often the 0.5 correction hides a real harm (true OR = 1.4)",
     cex.main = 0.95, font.main = 1)
axis(1, at = seq_len(nx), labels = paste0(ratios, ":1"), tick = FALSE)
axis(2, at = seq_len(ny), labels = paste0(rev(rarity), "x"), tick = FALSE)
for (i in seq_len(nrow(map))) {
  x <- match(map$ratio[i], ratios); y <- match(map$rarity[i], rev(rarity))
  z <- map$hidden_harm[i]
  col <- ramp[max(1, ceiling(100 * z / zmax))]
  rect(x - 0.48, y - 0.48, x + 0.48, y + 0.48, col = col, border = "white", lwd = 2)
  text(x, y, sprintf("%.0f%%", 100 * z), cex = 0.8,
       col = if (z / zmax > 0.55) "white" else ink)
}
points(pos_x(rosi$ratio), pos_y(rosi$rarity), pch = c(21, 24), bg = "#eb6834",
       col = "white", cex = 1.8, lwd = 1.5)
legend("bottom", legend = c("Rosiglitazone: heart attack", "Rosiglitazone: CV death"),
       pch = c(21, 24), pt.bg = "#eb6834", col = "white", pt.cex = 1.4, bty = "n",
       cex = 0.75, text.col = ink, horiz = TRUE, inset = c(0, -0.30), xpd = NA)
## colour key
par(mar = c(6.2, 0.5, 3, 3.2))
image(1, seq(0, zmax, length.out = 100), t(seq_len(100)), col = ramp, axes = FALSE,
      xlab = "", ylab = "")
axis(4, at = pretty(c(0, zmax)), labels = paste0(100 * pretty(c(0, zmax)), "%"),
     las = 1, cex.axis = 0.75, col.axis = ink)
mtext("Hidden harm", side = 3, line = 0.6, cex = 0.7, col = ink)
dev.off()

## ---- LaTeX table: fragility --------------------------------------------------
lab <- c(MH = "Mantel--Haenszel", MH05 = "MH + 0.5 correction", Peto = "Peto",
         UFirth = "Unconditional Firth", CFirth = "Conditional Firth")
cell <- function(o, m) { r <- frag[frag$outcome == o & frag$method == m, ]
  sprintf("%.3f & %s", r$p_original, ifelse(is.na(r$events_to_flip), "$>$50", r$events_to_flip)) }
rows <- sapply(names(lab), function(m)
  sprintf("%s & %s & %s \\\\", lab[m], cell("heart_attack", m), cell("cv_death", m)))
writeLines(c("\\begin{table}[t]",
  "\\caption{Fragility of each method's conclusion on the rosiglitazone data: the number of single-event changes needed to reverse statistical significance ($\\alpha = 0.05$).}",
  "\\label{tab:fragility}", "\\centering\\small", "\\begin{tabular}{lcccc}", "\\toprule",
  " & \\multicolumn{2}{c}{Heart attack} & \\multicolumn{2}{c}{Cardiovascular death} \\\\",
  "Method & $p$ & Events to flip & $p$ & Events to flip \\\\", "\\midrule", rows,
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "results/figures/table_fragility.tex")

cat("\nSaved results/map/map_results.csv, results/fragility.csv,",
    "results/figures/fig5_danger_map.pdf and results/figures/table_fragility.tex\n")
