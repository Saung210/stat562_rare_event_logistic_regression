# Reliable Inference in Rare-Event Meta-Analysis with Zero-Event Trials

Code for the STAT 562 (Fall 2026) project *Reliable Inference in Rare-Event
Meta-Analysis with Zero-Event Trials: A Comparative Study of Conditional Firth
Penalization and Standard Methods, with Application to Rosiglitazone*.

Authors: Samara Pires, Yuke Hou, Yanwen Su, Saung Hnin Phyu, Anna Tam Ly.

## What the code does

It compares five methods for estimating a common odds ratio from many trials
with rare events and frequent zero-event trials:

| Code name | Method |
|---|---|
| `MH`     | Mantel–Haenszel, no continuity correction (double-zero trials dropped) |
| `MH05`   | Mantel–Haenszel with 0.5 added to trials containing a zero |
| `Peto`   | Peto's one-step odds ratio |
| `UFirth` | Unconditional Firth logistic regression with trial intercepts (Evrenoglou et al., 2022) |
| `CFirth` | Conditional Firth: trial intercepts removed by conditioning, then Firth-penalized (this project) |

The methods are applied to the rosiglitazone meta-analysis (Nissen & Wolski,
2007; 42 trials) and evaluated in a Monte Carlo simulation that keeps the real
trial sizes and uses a known true odds ratio.

## Requirements

R (version 4.0 or later) and three packages:

```r
install.packages(c("metafor", "logistf", "metadat"))
```

`metadat` is optional: if it is missing (or older than version 1.4), the copy
of the data in `data/dat.nissen2007.rda` is used instead.

## Folder structure

```
README.md
data/dat.nissen2007.rda      rosiglitazone data (from the metadat package)
R/conditional_firth.R        conditional Firth estimator
R/methods.R                  the five methods, in a common format
R/simulate.R                 data generator, simulation loop, summary function
01_real_data.R               fit all methods to the real data
02_pilot.R                   pilot simulation (true OR = 1)
03_run_simulation.R          full simulation, all scenarios
04_summarise.R               bias, coverage, false-positive rate, power
05_figures.R                 figures and LaTeX table for the paper
06_danger_map.R              danger map, fragility index, rosiglitazone design
results/                     created automatically; all outputs go here
```

## How to run

Open R **in the project folder** (in RStudio: Session → Set Working Directory →
Choose Directory), then run the scripts in order:

```r
source("01_real_data.R")        # about 10 seconds
source("02_pilot.R")            # about a minute
source("03_run_simulation.R")   # roughly 20-60 minutes with R = 5000 (see below)
source("04_summarise.R")        # under a minute
source("05_figures.R")          # a few seconds
source("06_danger_map.R")       # roughly 10-20 minutes with R_MAP = 2000
```

`06_danger_map.R` (i) simulates a 5 x 5 grid of allocation ratio x event
rarity at true OR = 1.4 and records how often conditional Firth detects the
harm (p < 0.05) while the 0.5 correction does not ("hidden harm"); (ii)
reports the same rate at the real rosiglitazone design (from scenarios 03
and 04, so run it after `03_run_simulation.R`); (iii) computes a greedy
fragility index on the real data: the number of single-event changes needed
to reverse each method's statistical significance.

### Notes on the simulation

- **Test first.** Set `R <- 20` at the top of `03_run_simulation.R`, run it,
  delete the `results/sim` folder, then set `R <- 5000` for the real run.
- **Run time.** About 0.06 seconds per replicate per core. With 8 cores,
  all 14 scenarios at R = 5000 take roughly 20 minutes; with 2 cores about
  an hour.
- **Parallel.** Scripts use all cores but two (`parallel::detectCores() - 2`),
  on Windows, Mac and Linux. Results are identical whatever the number of
  cores, because every replicate has its own fixed seed.
- **Restarting.** Each scenario is saved with its settings as soon as it
  finishes. On a restart, a scenario is skipped only if its saved settings
  (including R and the seed) match; otherwise it is rerun, so test runs are
  never mixed with the real run.
- **Baseline risk.** When control risk varies between trials (`p0_sd = 1`),
  the logit-normal mean is calibrated so the average control risk stays equal
  to the pooled rate (`calibrate_mu()` in `R/simulate.R`).
- **Failed fits.** A fit counts as successful only if the estimate, both CI
  limits and the p-value are finite. Rates are computed among successful fits
  (`n_ok`); `n_failed` and `n_warning` are reported alongside. The warnings
  recorded for MH and Peto are metafor's notice that zero-cell trials were
  omitted from the trial-level estimates, and are expected.
- **Reproducibility.** Scenario seeds are fixed in `03_run_simulation.R`; the
  scenario settings are written to `results/scenarios.csv`.

## Outputs

| File | Contents |
|---|---|
| `results/real_data_summary.csv` | zero-event trials, event rates, allocation |
| `results/real_data_results.csv` | OR, 95% CI and p-value from every method |
| `results/pilot_summary.csv` | pilot false-positive rate and coverage |
| `results/sim/scenario_XX.rds` | raw simulation output, one file per scenario |
| `results/summary_all.csv` | bias, coverage, rejection rate, Monte Carlo SEs |
| `results/firth_comparison.csv` | paired unconditional vs conditional Firth |
| `results/figures/` | Figures 1–5 (PDF) and LaTeX tables (real data, fragility) |
| `results/map/map_results.csv` | hidden-harm rate, power and bias for each map cell |
| `results/map/real_design_hidden_harm.csv` | hidden-harm rate at the real design |
| `results/fragility.csv` | fragility index for each method and outcome |

All simulation summaries report Monte Carlo standard errors.

## Data source

Nissen, S. E., & Wolski, K. (2007). Effect of rosiglitazone on the risk of
myocardial infarction and death from cardiovascular causes. *New England
Journal of Medicine, 356*(24), 2457–2471. Data as distributed in the R package
`metadat` (`dat.nissen2007`).
