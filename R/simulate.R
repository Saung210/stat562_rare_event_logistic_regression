
# R/simulate.R
# Data generation, simulation, and summary functions.
# Run R/methods.R first.

library(parallel)

# Calibrate the mean control risk
calibrate_mu <- function(p0, p0_sd) {
  if (p0_sd == 0) return(qlogis(p0))
  
  mean_risk <- function(mu) {
    integrate(function(z) plogis(mu + p0_sd * z) * dnorm(z),
              -Inf, Inf, rel.tol = 1e-10)$value
  }
  
  uniroot(function(mu) mean_risk(mu) - p0,
          interval = qlogis(p0) + c(-10, 2), tol = 1e-12)$root
}

# Generate data for one simulated meta-analysis
gen_data <- function(n1, n0, p0, OR, p0_sd = 0, mu = NULL) {
  K <- length(n1)
  
  if (p0_sd > 0) {
    if (is.null(mu)) mu <- calibrate_mu(p0, p0_sd)
    p0 <- plogis(mu + rnorm(K, 0, p0_sd))
  } else {
    p0 <- rep_len(p0, K)
  }
  
  # Calculate treatment risk using the true OR
  p1 <- plogis(qlogis(p0) + log(OR))
  
  list(
    x1 = rbinom(K, n1, p1),
    n1 = n1,
    x0 = rbinom(K, n0, p0),
    n0 = n0
  )
}

# Fit all selected methods
fit_all <- function(s, method_names = names(methods)) {
  do.call(rbind, lapply(method_names, function(nm) {
    msgs <- character(0)
    err <- NA_character_
    
    r <- tryCatch(
      withCallingHandlers(
        methods[[nm]](s$x1, s$n1, s$x0, s$n0),
        warning = function(w) {
          msgs <<- c(msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        err <<- conditionMessage(e)
        c(logOR = NA, lo = NA, hi = NA, p = NA)
      }
    )
    
    # Check whether the fit returned valid results
    valid <- all(is.finite(r[c("logOR", "lo", "hi", "p")]))
    
    data.frame(
      method = nm,
      logOR = r[["logOR"]],
      lo = r[["lo"]],
      hi = r[["hi"]],
      p = r[["p"]],
      error = !valid,
      warning = length(msgs) > 0,
      warning_text = if (length(msgs)) {
        paste(unique(msgs), collapse = " | ")
      } else NA,
      error_text = err
    )
  }))
}

# Run simulations for one scenario
run_sim <- function(R, OR, p0, n1, n0, seed, cores = 1, p0_sd = 0,
                    method_names = names(methods)) {
  
  force(OR)
  force(p0)
  force(n1)
  force(n0)
  force(seed)
  force(p0_sd)
  force(method_names)
  
  # Calculate mu once for each scenario
  mu <- calibrate_mu(p0, p0_sd)
  
  one_rep <- function(r) {
    set.seed(seed * 100000 + r)
    
    s <- gen_data(n1, n0, p0, OR, p0_sd, mu)
    res <- fit_all(s, method_names)
    
    res$rep <- r
    res$zero_any <- sum(s$x1 == 0 | s$x0 == 0)
    res$zero_double <- sum(s$x1 == 0 & s$x0 == 0)
    res$events <- sum(s$x1 + s$x0)
    
    res
  }
  
  # Use parallel processing when more than one core is available
  if (cores <= 1) {
    out <- lapply(seq_len(R), one_rep)
    
  } else if (.Platform$OS.type == "windows") {
    cl <- makeCluster(cores)
    on.exit(stopCluster(cl), add = TRUE)
    
    clusterCall(cl, setwd, getwd())
    clusterEvalQ(cl, {
      suppressMessages(source("R/methods.R"))
      source("R/simulate.R")
    })
    
    out <- parLapply(cl, seq_len(R), one_rep)
    
  } else {
    out <- mclapply(seq_len(R), one_rep, mc.cores = cores)
  }
  
  out <- do.call(rbind, out)
  out$true_logOR <- log(OR)
  
  out
}

# Summarize the simulation results
summarise_sim <- function(sim, alpha = 0.05) {
  do.call(rbind, lapply(split(sim, sim$method), function(s) {
    truth <- s$true_logOR[1]
    ok <- !s$error
    est <- s$logOR[ok]
    n <- sum(ok)
    
    cover <- if (n > 0) {
      mean(s$lo[ok] <= truth & truth <= s$hi[ok])
    } else NA
    
    reject <- if (n > 0) {
      mean(s$p[ok] < alpha)
    } else NA
    
    data.frame(
      method = s$method[1],
      n_ok = n,
      n_failed = sum(!ok),
      n_warning = sum(s$warning),
      bias = if (n > 0) mean(est) - truth else NA,
      bias_mcse = if (n > 1) sd(est) / sqrt(n) else NA,
      emp_se = if (n > 1) sd(est) else NA,
      coverage = cover,
      cover_mcse = sqrt(cover * (1 - cover) / n),
      reject = reject,
      reject_mcse = sqrt(reject * (1 - reject) / n)
    )
  }))
}
