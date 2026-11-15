
# R/methods.R
# Five methods for estimating the odds ratio.
# Each function returns the log OR, 95% CI, and two-sided p-value.

library(metafor)
library(logistf)
source("R/conditional_firth.R")

out <- function(est, lo, hi, p) {
  c(logOR = unname(est), lo = unname(lo), hi = unname(hi), p = unname(p))
}

# 1. Mantel-Haenszel without continuity correction
m_mh <- function(x1, n1, x0, n0) {
  m <- rma.mh(ai = x1, n1i = n1, ci = x0, n2i = n0, measure = "OR",
              add = c(0.5, 0), to = c("only0", "none"),
              drop00 = c(TRUE, TRUE))
  out(m$beta[1], m$ci.lb, m$ci.ub, m$pval)
}

# 2. Mantel-Haenszel with 0.5 continuity correction
m_mh05 <- function(x1, n1, x0, n0) {
  m <- rma.mh(ai = x1, n1i = n1, ci = x0, n2i = n0, measure = "OR",
              add = c(0.5, 0.5), to = c("only0", "only0"),
              drop00 = c(TRUE, TRUE))
  out(m$beta[1], m$ci.lb, m$ci.ub, m$pval)
}

# 3. Peto's method
m_peto <- function(x1, n1, x0, n0) {
  m <- rma.peto(ai = x1, n1i = n1, ci = x0, n2i = n0)
  out(m$beta[1], m$ci.lb, m$ci.ub, m$pval)
}

# 4. Unconditional Firth logistic regression
m_ufirth <- function(x1, n1, x0, n0) {
  K <- length(x1)
  
  # Set up the data for logistic regression
  long <- data.frame(
    trial = factor(rep(seq_len(K), 4)),
    treat = rep(c(1, 1, 0, 0), each = K),
    y = rep(c(1, 0, 1, 0), each = K),
    w = c(x1, n1 - x1, x0, n0 - x0)
  )
  
  long <- long[long$w > 0, ]
  long$trial <- droplevels(long$trial)
  
  # Increase the iteration limit and calculate the CI for treatment only
  ctrl <- logistf.control(maxit = 100)
  pctrl <- logistpl.control(maxit = 100)
  
  m <- withCallingHandlers(
    logistf(y ~ treat + trial, data = long, weights = w,
            plconf = 2, control = ctrl, plcontrol = pctrl),
    warning = function(w) {
      if (grepl("iteration|converg", conditionMessage(w), ignore.case = TRUE))
        stop("logistf did not converge: ", conditionMessage(w))
    }
  )
  
  stopifnot(names(coef(m))[2] == "treat")
  out(coef(m)[2], m$ci.lower[2], m$ci.upper[2], m$prob[2])
}

# 5. Conditional Firth method
m_cfirth <- function(x1, n1, x0, n0) {
  m <- cond_firth(x1, n1, x0, n0)
  
  # Treat estimates at the search boundary as infinite
  est <- if (m$at.bound) sign(m$logOR) * Inf else m$logOR
  
  out(est, log(m$ci.OR[1]), log(m$ci.OR[2]), m$p.value)
}

# List of methods used in the simulation
methods <- list(
  MH = m_mh,
  MH05 = m_mh05,
  Peto = m_peto,
  UFirth = m_ufirth,
  CFirth = m_cfirth
)
