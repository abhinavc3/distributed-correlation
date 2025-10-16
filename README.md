
# distributed-correlation
Collects the code for the simulation section of the paper *“When Data Can’t Meet: Estimating Correlation Across Privacy Barriers.”*

## Repository layout
- **`vert-cor.R`** — Main simulation driver for differentially private correlation estimators (non-interactive and interactive variants). Sets seeds, grids over distributions (e.g., Gaussian/Bernoulli), and runs Monte Carlo replications. Also contains small utility routines (e.g., mixture quantiles) and supports parallelism via `parallel::mclapply`.  
- **`ver-cor-subG.R`** — Core methods for **sub-Gaussian-calibrated** DP correlation inference. Implements confidence-interval constructors (NI and INT) with clipping/threshold choices (`lambda_n`, `lambda_INT_n`) and helper routines (e.g., DP mean / DP second-moment pieces).  
- **`real-data-sims.R`** — Reproducible **HRS** real-data experiment (BMI vs. Age) under DP. Loads the supplied RDS panel, summarizes missingness by wave, clips, and computes private means/variances and DP correlation; also includes minimalist plotting with `ggplot2`.  
- **`hrs_long_panel.rds`** — Prepared HRS long-format panel used by `real-data-sims.R` (wave-level observations with BMI / age and metadata).  

> Note: Files prefixed by `__MACOSX/` are macOS archive artifacts and can be ignored.

## Quick start
**R version:** R ≥ 4.1 recommended.

**Required packages:**
- `extraDistr` (Laplace noise helpers)  
- `dplyr`, `tidyr`, `ggplot2` (real-data script)  
- `parallel` (optional, for faster sims on Unix/Mac)

To install:
```r
install.packages(c("extraDistr", "dplyr", "tidyr", "ggplot2"))


## How to run

### Simulations (synthetic)

Edit grids/seeds in `vert-cor.R` if desired, then run:

```r
source("vert-cor.R")
```

Results (estimates / CIs) are produced in-session; adapt the code to save tables or plots if needed.

### HRS real data

Ensure `hrs_long_panel.rds` is in the working directory and run:

```r
source("real-data-sims.R")
```

This computes DP summaries and produces simple diagnostic plots.

## Notes

* Privacy parameters (ε’s) and clipping ranges are set inside each script; adjust them to match the paper’s settings.
* The sub-Gaussian tuning (`lambda_n`, `lambda_INT_n`) in `ver-cor-subG.R` mirrors the paper’s theory; comments in the file mark any deliberate deviations for robustness.

```
```
