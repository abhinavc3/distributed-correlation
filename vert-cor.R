############################################################
##  Simulation study: Private correlation estimation      ##
##  Author: <your name>                                   ##
##  Date:   2025-05-04                                    ##
############################################################

## ---- 0.1  Packages --------------------------------------
# All you need is base R + extraDistr for Laplace noise.
if (!requireNamespace("extraDistr", quietly = TRUE))
  install.packages("extraDistr")

library(extraDistr)  # dlaplace / rlaplace helpers
library(parallel)    # for mclapply on *nix

## ---- 0.2  Global reproducibility ------------------------
MASTER_SEED <- 2025L
set.seed(MASTER_SEED)

## ---- 0.3  Simulation grid -------------------------------
# Sample sizes
n_grid   <- c(200, 400, 800, 1600, 3200)

# Privacy-budget pairs (ε₁ serverX, ε₂ serverY)
eps_grid <- rbind(
  c(0.2, 0.2),
  c(0.5, 0.5),
  c(1.0, 1.0),
  c(1.5, 0.5),
  c(0.5, 1.5)
)
colnames(eps_grid) <- c("eps1", "eps2")

# Correlation coefficients
rho_grid <- c(0, 0.3, 0.8)

# Distributions to simulate
distr_grid <- c("gaussian", "bernoulli")

# Monte-Carlo replications per design point
B <- 1000

## ---- 0.4  Quantiles of Mixture of Gaussian and scaled Laplace------

mixquant<- function(c,p){
  #set.seed(1)
  nsim <- 1000
  xvec <- rnorm(nsim) + c*rexp(nsim)*(2*rbinom(nsim,1,0.5)-1)
  return(sort(xvec)[ceiling(p*nsim)])
  
  #xgrid <- seq(-10,10,0.01)
  #ygrid <- seq(-10,10,0.01)
  #prod <- function(x){dnorm(x-ygrid)*exp(-abs(ygrid)/c)/(2*c)*0.01} 
  #pdfmix <- rowSums(sapply(xgrid, prod))
  #cdfmix <- cumsum(pdfmix)*0.01
  #return(min(xgrid[cdfmix>p]))
}


## =========================================================
## 1. DATA-GENERATION FUNCTIONS                            ##
## =========================================================

# 1.1  Gaussian with target correlation rho ---------------
gen_gaussian <- function(n, rho, mu = c(0, 0)) {
  ## Returns an n×2 matrix with columns (X, Y).
  Sigma <- matrix(c(1, rho,  # var(X)=1, cov(X,Y)=rho
                    rho, 1), # symmetric
                  nrow = 2, byrow = TRUE)
  # Use built-in MASS::mvrnorm to avoid writing our own.  
  if (!requireNamespace("MASS", quietly = TRUE))
    install.packages("MASS")
  MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
}

# 1.2  Correlated Bernoulli (marginals 0.5) ---------------
#  - X,Y ∈ {0,1} with P(X=1)=P(Y=1)=0.5 and Corr(X,Y)=rho
#  - Valid for rho ∈ [-1, 1].
gen_bernoulli <- function(n, rho) {
  stopifnot(abs(rho) <= 1)
  p11 <- 0.25 + rho / 4  # see derivation in notes
  p10 <- 0.25 - rho / 4
  p01 <- p10
  p00 <- p11            # symmetry ensures p00 = p11
  ## CDF inversion for joint sample -----------------------
  u <- runif(n)
  v <- runif(n)
  # Draw X first
  X <- as.numeric(u < 0.5)
  # Conditional law of Y given X
  Y <- numeric(n)
  # Case X=0
  idx0 <- which(X == 0)
  Y[idx0] <- as.numeric(v[idx0] < (p01 / 0.5))
  # Case X=1
  idx1 <- which(X == 1)
  Y[idx1] <- as.numeric(v[idx1] < (p11 / 0.5))
  cbind(X, Y)
}

## =========================================================
## 2.  SIGN–BATCH PRIVATE ESTIMATOR  (Gaussian only)       ##
##     Source: Sec. 3.1, eqs. (1)–(2) of the paper         ##
## =========================================================

# Laplace generator --------------------------------------------------------
rLap <- function(n = 1, scale) extraDistr::rlaplace(n, mu = 0, sigma = scale)

# 2.1  Non-interactive sign–batch estimator -------------------------------
#
# Args:
#   X, Y   : numeric vectors of equal length n (Gaussian, mean 0, var 1)
#   eps1   : DP budget for server X
#   eps2   : DP budget for server Y
# Returns: scalar ρ̂ (private correlation estimate)
#

# X and Y should be noramlised before passing to this function
correlation_NI_signbatch <- function(X, Y, eps1, eps2) {
  stopifnot(length(X) == length(Y), eps1 > 0, eps2 > 0)
  n <- length(X)
  

  ## ----------  Step 1: choose batch size  m  ----------------------------
  m  <- ceiling(8 / (eps1 * eps2))          # paper’s optimal choice
  if (m > n) m <- n                         # fallback if n is tiny
  k  <- floor(n / m)                        # number of full batches
  if (k < 1) stop("Need at least one full batch; increase n or reduce m")
  
  ## ----------  Step 2: aggregate signed data in each batch --------------
  # (use only the first  k*m  observations to keep code simple)
  idx_use <- seq_len(k * m)
  Xs <- sign(X[idx_use])
  Ys <- sign(Y[idx_use])
  
  X_bar <- Y_bar <- numeric(k)
  for (j in seq_len(k)) {
    idx <- ((j - 1) * m + 1):(j * m)
    X_bar[j] <- mean(Xs[idx])
    Y_bar[j] <- mean(Ys[idx])
  }
  
  ## ----------  Step 3: add Laplace noise to each batch mean -------------
  # Sensitivity of a batch mean of signs is  2/m,
  # so Laplace scale parameter = (2/m) / ε   (eq. (1))
  X_noisy <- X_bar + rLap(k, scale = 2 / (m * eps1))
  Y_noisy <- Y_bar + rLap(k, scale = 2 / (m * eps2))
  
  ## ----------  Step 4: form η̂ and back-transform to ρ̂  -----------------
  # eq. (2):   η̂ = (m / k) * Σ_j  X̃_j Ỹ_j
  eta_hat <- (m / k) * sum(X_noisy * Y_noisy)
  
  # eq. (2) continued:  ρ̂ = sin(π η̂ / 2)
  rho_hat <- sin(pi * eta_hat / 2)
  
  return(rho_hat)
}


## =========================================================
## 2A.  GENERIC ONE-STEP INTERACTIVE ESTIMATOR (Gaussian)  ##
##       – chooses the “sender’’ by the larger ε           ##
## =========================================================
# X and Y should be noramlised before passing to this function
correlation_INT_signflip <- function(X, Y, eps1, eps2) {
  stopifnot(length(X) == length(Y), eps1 > 0, eps2 > 0)
  n <- length(X)
  

  ## ---- 1. Decide roles ----------------------------------
  sender_is_X <- (eps1 >= eps2)        # TRUE  ⇢  X → Y   (ε_s = ε1)
  eps_s <- if (sender_is_X) eps1 else eps2   # larger ε  (sender)
  eps_r <- if (sender_is_X) eps2 else eps1   # smaller ε (receiver)
  
  p <- exp(eps_s) / (exp(eps_s) + 1)          # flip prob
  S <- rbinom(n, 1, p)                        # Bernoulli flips (1 keeps sign)
  
  ## ---- 2. Build the privacy-protected sum ---------------
  if (sender_is_X) {
    core <- (2 * S - 1) * sign(X) * sign(Y)   # flip sign(X)
  } else {
    core <- (2 * S - 1) * sign(Y) * sign(X)   # flip sign(Y)
  }
  sum_core <- sum(core)
  
  ## ---- 3. Add Laplace noise for receiver privacy --------
  scale_Z <- 2 * (exp(eps_s) + 1) /
    (n * (exp(eps_s) - 1) * eps_r)
  Z <- extraDistr::rlaplace(1, 0, scale_Z)
  
  eta_hat <- (exp(eps_s) + 1) /
    (n * (exp(eps_s) - 1)) * sum_core + Z
  
  ## ---- 4. Back-transform to ρ̂ ---------------------------
  sin(pi * eta_hat / 2)
}

## =========================================================
## 2B.  CONFIDENCE INTERVALS                                ##
##      – NI sign-batch (Sec 3.1)                           ##
##      – INT sign-flip  (Sec 4.1.1)                        ##
## =========================================================

## 2B.1  NI sign-batch CI  ---------------------------------
ci_NI_signbatch <- function(X, Y, eps1, eps2, alpha = 0.05, normalise = T) {
  
  n <- length(X)
  m <- ceiling(8 / (eps1 * eps2))
  k <- floor(n / m)
  stopifnot(k >= 1)
  
  if(normalise == T){
    L_clip = sqrt(2*log(n))
    ##  private normalise + clip (RAW Gaussians, *not* the signs)
    X <- priv_standardize(X, eps1, L_clip)
    Y <- priv_standardize(Y, eps2, L_clip)
  }
  idx_use <- seq_len(k * m)
  Xs <- sign(X[idx_use])
  Ys <- sign(Y[idx_use])
  
  # batch means of signs
  X_bar <- Y_bar <- numeric(k)
  for (j in seq_len(k)) {
    idx <- ((j - 1) * m + 1):(j * m)
    X_bar[j] <- mean(Xs[idx])
    Y_bar[j] <- mean(Ys[idx])
  }
  
  # Laplace-noised batch means
  X_tilde <- X_bar + rLap(k, scale = 2 / (m * eps1))
  Y_tilde <- Y_bar + rLap(k, scale = 2 / (m * eps2))
  
  Tj <-    m* X_tilde * Y_tilde                    # Sec 3.1, eq. (2) component
  eta_hat <- (1/ k) * sum(Tj)
  rho_hat <- sin(pi * eta_hat / 2)
  
  ## asymptotic SE  (Sec 3.1 Asymp. dist.) -----------------
  
  S_eta <- sd(Tj)
  #S_eta <- sqrt(((m-1)/m)^2*(1+eta_hat^2)+1/m+8*(1/eps1^2+1/eps2^2)/m+64/(m*eps1*eps2)^2)  # Plug in estimate of E[eta]
  se <- (pi / 2) * S_eta * sqrt(1 - rho_hat^2) / sqrt(k)
  crit <- qnorm(1 - alpha / 2)
  
  #list(
  #  rho_hat = rho_hat,
  #  ci      = c(max(rho_hat - crit * se,-1),
  #              min(rho_hat + crit * se, 1))
  #)
  list(
    rho_hat = rho_hat,
    ci = c(sin(pi/2*max( eta_hat - crit * S_eta / sqrt(k), -1)),
           sin(pi/2*min(eta_hat + crit * S_eta / sqrt(k), 1)) 
           )
  )
}

## =========================================================
## 2B.2  CI for INTERACTIVE SIGN-FLIP (role-agnostic)      ##
## =========================================================
ci_INT_signflip <- function(X, Y, eps1, eps2,
                            alpha = 0.05,
                            mode = c("auto", "normal", "laplace"),
                            normalise = T) {
  stopifnot(length(X) == length(Y), eps1 > 0, eps2 > 0)
  n <- length(X)
  mode <- match.arg(mode)
  
  if(normalise == T){
    L_clip = sqrt(2*log(n))
    ##  private normalise + clip (RAW Gaussians, *not* the signs)
    X <- priv_standardize(X, eps1, L_clip)
    Y <- priv_standardize(Y, eps2, L_clip)
  }
  ## ---- 1. Role assignment (same as estimator) -----------
  sender_is_X <- (eps1 >= eps2)
  eps_s <- if (sender_is_X) eps1 else eps2   # sender ε
  eps_r <- if (sender_is_X) eps2 else eps1   # receiver ε
  
  ## ---- 2. Point estimate --------------------------------
  rho_hat <- correlation_INT_signflip(X, Y, eps1, eps2)
  eta_hat <- 1-acos(rho_hat)*2/pi
  
  ## ---- 3. Variance of η̂ (depends only on eps_s, eps_r) --
  sigma_eta2 <- 1- ((exp(eps_s)-1)/(exp(eps_s)+1))^2*(1-acos(rho_hat)*2/pi)^2
    #4 * exp(eps_s) / (exp(eps_s) + 1)^2
  #scale_Z    <- 2 * (exp(eps_s) + 1) /
  #  (n * (exp(eps_s) - 1) * eps_r)
  #var_eta    <- sigma_eta2 / n #+ 2 * scale_Z^2   # Var Laplace = 2·scale²
  ratio <- (exp(eps_s) + 1) / (exp(eps_s) - 1)
  se_norm    <- 1/sqrt(n) * (pi / 2) * sqrt(1 - rho_hat^2) * sqrt(sigma_eta2) * ratio
  se_norm_eta <- 1/sqrt(n)* sqrt(sigma_eta2) * ratio
  
  ## ---- 4. Normal or Laplace dominated interval ----------
  if (mode == "auto") {
    mode <- if (sqrt(n) * eps_r > 0.5) "normal" else "laplace"
  }
  
  if (mode == "normal") {            # Case 1 in §4.1.1
    cstar <- 2/(sqrt(n*sigma_eta2)*eps_r)
    #width <- mixquant(cstar, 1-alpha/2) * se_norm
      #qnorm(1 - alpha / 2) * se_norm
    width_eta <- mixquant(cstar, 1-alpha/2) * se_norm_eta
  } else {                           # Case 2 in §4.1.1
    ratio <- (exp(eps_s) + 1) / (exp(eps_s) - 1)
    scale_L <- (pi * sqrt(1 - rho_hat^2) / (n * eps_r)) * ratio
    scale_L_eta <- (2 / (n * eps_r)) * ratio
    #width <- scale_L * log(1 / alpha)
    width_eta <- scale_L_eta * log(1 / alpha)
  }
  
  list(rho_hat = rho_hat,
       ci      = c(sin(pi/2*max(eta_hat - width_eta, -1)),
                   sin(pi/2*min(eta_hat + width_eta, 1))),
       #ci      = c(rho_hat - width_eta, rho_hat + width),
       mode    = mode,
       roles   = if (sender_is_X) "X→Y" else "Y→X")
}

############################################################################
##  PRIVATE  centre–scale  (single pre-clip, no post-clip)                ##
############################################################################
priv_standardize <- function(vec,
                             eps_norm,
                             L_raw = 6) {   # clip once, before adding noise
  n <- length(vec)
  
  #### 1) hard-clip raw values  -------------------------------------------
  x_clipped <- pmax(pmin(vec,  L_raw), -L_raw)
  
  #### 2) split ε equally for μ and M₂  -----------------------------------
  eps_mu <- eps_norm / 2
  eps_m2 <- eps_norm / 2
  
  ## DP mean  — sensitivity = 2 L_raw / n
  mu_priv <- mean(x_clipped) +
    rLap(1, scale = 2 * L_raw / (n * eps_mu))
  
  ## DP second moment M₂ = mean(X²) — sensitivity = 2 L_raw² / n
  m2_priv <- mean(x_clipped^2) +
    rLap(1, scale = 2 * L_raw^2 / (n * eps_m2))
  
  ## DP variance & SD  (ensure non-negative)
  var_priv <- max(m2_priv - mu_priv^2, 1e-12)
  sd_priv  <- sqrt(var_priv)
  
  #### 3) standardise (no further clipping) -------------------------------
  (x_clipped - mu_priv) / sd_priv
}


## =========================================================
## 4.  PARAMETER–SWEEP DRIVER                              ##
##     – Gaussian only, NI sign–batch + INT sign-flip      ##
## =========================================================

run_sim_one <- function(n, rho, eps1, eps2, 
                        mu = c(0, 0), sigma = c(1, 1),
                        B      = 1000,
                        alpha  = 0.05,
                        ci_mode = "auto",
                        normalise = T,
                        seed   = 2025L) {
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
  set.seed(seed)                       # reproducibility
  
  ## ---- storage ------------------------------------------
  out <- data.frame(
    repl        = seq_len(B),
    # point estimates
    ni_hat      = NA_real_,
    int_hat     = NA_real_,
    # squared errors
    ni_se2      = NA_real_,
    int_se2     = NA_real_,
    # CI bounds  (NEW)
    ni_low      = NA_real_,
    ni_up       = NA_real_,
    int_low     = NA_real_,
    int_up      = NA_real_,
    # convenience: coverage flags & length
    ni_cover    = NA_integer_,
    int_cover   = NA_integer_,
    ni_ci_len   = NA_real_,
    int_ci_len  = NA_real_
  )
  
  
  ## ---- Monte-Carlo loop ---------------------------------
  Sigma <- matrix(c(sigma[1]^2, sigma[1]*sigma[2]*rho,
                    sigma[1]*sigma[2]*rho, sigma[2]^2), 2, 2)
  
  for (b in seq_len(B)) {
    ## 1.  Generate Gaussian sample ------------------------
    XY <- MASS::mvrnorm(n, mu, Sigma = Sigma)
    X  <- XY[, 1];  Y <- XY[, 2]
    
    ## ---- 2.  Non-interactive estimate + CI -------------------
    ni_res <- ci_NI_signbatch(X, Y, eps1, eps2,
                              alpha = alpha, normalise = normalise)
    
    out$ni_hat [b] <- ni_res$rho_hat
    out$ni_low [b] <- ni_res$ci[1]
    out$ni_up  [b] <- ni_res$ci[2]
    out$ni_se2 [b] <- (ni_res$rho_hat - rho)^2
    out$ni_cover[b] <- rho >= ni_res$ci[1] && rho <= ni_res$ci[2]
    out$ni_ci_len[b] <- diff(ni_res$ci)
    
    ## ---- 3.  Interactive estimate + CI -----------------------
    int_res <- ci_INT_signflip(X, Y, eps1, eps2,
                               alpha = alpha, mode = ci_mode, normalise = normalise)
    
    out$int_hat [b] <- int_res$rho_hat
    out$int_low [b] <- int_res$ci[1]
    out$int_up  [b] <- int_res$ci[2]
    out$int_se2 [b] <- (int_res$rho_hat - rho)^2
    out$int_cover[b] <- rho >= int_res$ci[1] && rho <= int_res$ci[2]
    out$int_ci_len[b] <- diff(int_res$ci)
    
  }
  
  ## ---- aggregate summaries ------------------------------
  summarise <- function(est, se2, cover, lo, up) {
    c(
      mse       = mean(se2),
      bias      = mean(est) - rho,
      var       = var(est),
      coverage  = mean(cover),
      ci_length = mean(up - lo)
    )
  }
  
  summ_df <- rbind(
    NI  = summarise(out$ni_hat,  out$ni_se2,
                    out$ni_cover, out$ni_low,  out$ni_up),
    INT = summarise(out$int_hat, out$int_se2,
                    out$int_cover, out$int_low, out$int_up)
  )
  
  summ_df <- as.data.frame(summ_df)
  summ_df$method <- rownames(summ_df)
  rownames(summ_df) <- NULL
  
  list(detail = out, summary = summ_df)
}



## single design point ------------------------------------
res <- run_sim_one(
  n     = 2000,
  rho   = -0.95,
  eps1  = 0.5,
  eps2  = 1,
  mu = c(2, 2),
  sigma = c(2, 0.1),
  normalise = T,
  B     = 1000        # increase for final runs
)

print(res$summary)
#>   mse        bias        var coverage ci_length method
#> 1 ...  (non-interactive stats)
#> 2 ...  (interactive   stats)

## full replicate matrix is in res$detail
head(res$detail)










################################################################
##  5.  GRID-WISE SIMULATION + FIGURES                        ##
##      – Gaussian only, NI sign–batch & INT sign-flip        ##
##      – Requires run_sim_one() from Chunk 4                 ##
################################################################

## ------------------------------------------------------------
## 5.1  Define the design grid                                ##
## ------------------------------------------------------------
mu <- c(0.5,0.5)
sigma <- c(2,2)
n_grid   <- c(1000, 1500, 2500, 4000, 6000, 9000)
rho_grid <- c(0, 0.15, 0.3, 0.4, 0.5, 0.65, 0.8, 0.9)
eps_pairs <- list(
  c(0.5, 0.5),
  c(1.0, 1.0),
  c(1.5, 0.5)     # asymmetric
)

B_per_setting <- 250               # ↑ for camera-ready
alpha_level   <- 0.05
ci_mode_choice <- "auto"
normalise <- T
## ------------------------------------------------------------
## 5.2  Run the whole grid  (~ minutes not hours)             ##
## ------------------------------------------------------------
expand_grid <- function(...) {
  expand.grid(..., KEEP.OUT.ATTRS = FALSE)
}

design_df <- expand_grid(
  n      = n_grid,
  rho    = rho_grid,
  eps_idx= seq_along(eps_pairs)      # just an index into eps_pairs
)

n_cores <- max(1L, parallel::detectCores() - 1L)

if (.Platform$OS.type == "windows") {
  message("Windows detected – running serially (mc.cores set to 1).")
  grid_results <- lapply(seq_len(nrow(design_df)), function(i) {
    row <- design_df[i, ]
    eps1 <- eps_pairs[[row$eps_idx]][1]
    eps2 <- eps_pairs[[row$eps_idx]][2]
    
    run_sim_one(n      = row$n,
                rho    = row$rho,
                eps1   = eps1,
                eps2   = eps2,
                mu     = mu,
                sigma  = sigma,
                B      = B_per_setting,
                alpha  = alpha_level,
                ci_mode= ci_mode_choice,
                seed   = 1e6 + i)
  })
} else {                               # Linux / macOS – keep multicore
  grid_results <- parallel::mclapply(
    seq_len(nrow(design_df)),
    mc.cores = n_cores,
    FUN = function(i) {
      row <- design_df[i, ]
      eps1 <- eps_pairs[[row$eps_idx]][1]
      eps2 <- eps_pairs[[row$eps_idx]][2]
      
      run_sim_one(n      = row$n,
                  rho    = row$rho,
                  eps1   = eps1,
                  eps2   = eps2,
                  mu     = mu,
                  sigma  = sigma,
                  normalise = normalise,
                  B      = B_per_setting,
                  alpha  = alpha_level,
                  ci_mode= ci_mode_choice,
                  seed   = 1e6 + i)
    })
}

## 5.2.1  Merge replicate-level data --------------------------
for (i in seq_len(nrow(design_df))) {
  df_i <- design_df[i, ]
  detail <- grid_results[[i]]$detail
  detail$n        <- df_i$n
  detail$rho_true <- df_i$rho
  detail$eps1     <- eps_pairs[[df_i$eps_idx]][1]
  detail$eps2     <- eps_pairs[[df_i$eps_idx]][2]
  grid_results[[i]]$detail <- detail
}

library(data.table)
detail_all <- rbindlist(lapply(grid_results, `[[`, "detail"))
saveRDS(detail_all, "sim_detail_all.rds")   # raw data backup


library(data.table)

## ----- INT summary ---------------------------------------
summ_INT <- detail_all[
  , .(
    mse       = mean(int_se2),
    bias      = mean(int_hat) - mean(rho_true),
    coverage  = mean(int_cover),
    ci_len    = mean(int_ci_len)
  ),
  by = .(n, rho_true, eps1, eps2)
][, method := "INT"]      # add a constant column afterwards

## ----- NI summary ----------------------------------------
summ_NI <- detail_all[
  , .(
    mse       = mean(ni_se2),
    bias      = mean(ni_hat) - mean(rho_true),
    coverage  = mean(ni_cover),
    ci_len    = mean(ni_ci_len)
  ),
  by = .(n, rho_true, eps1, eps2)
][, method := "NI"]

## combine -------------------------------------------------
summ_all <- rbindlist(list(summ_NI, summ_INT))


## ------------------------------------------------------------
## 5.3  Figure 1 – Mean offset bands vs ρ  (centre = ρ̂ − ρ)  ##
## ------------------------------------------------------------
library(data.table)
library(ggplot2)

chosen_n   <- 1500
chosen_eps <- c(1.5, 0.5)

slice <- detail_all[
  n   == chosen_n &
    eps1 == chosen_eps[1] &
    eps2 == chosen_eps[2]
]

## --- mean lower / upper offsets *and* mean ρ̂ offset -------
band_dt <- rbind(
  slice[, .(
    rho        = rho_true,
    lower_off  = mean(ni_low  - rho_true),
    upper_off  = mean(ni_up   - rho_true),
    est_off    = mean(ni_hat  - rho_true),   # <- NEW
    method     = "NI"
  ), by = rho_true],
  slice[, .(
    rho        = rho_true,
    lower_off  = mean(int_low - rho_true),
    upper_off  = mean(int_up  - rho_true),
    est_off    = mean(int_hat - rho_true),   # <- NEW
    method     = "INT"
  ), by = rho_true]
)

band_dt[, method := factor(method, levels = c("NI","INT"))]

title_expr <- bquote(
  "Mean CI offset bands - n = " * .(chosen_n) *
    ",  " * epsilon[1] * " = " * .(chosen_eps[1]) *
    ",  " * epsilon[2] * " = " * .(chosen_eps[2])
)

p_band <- ggplot(band_dt, aes(x = rho, fill = method)) +
  geom_ribbon(aes(ymin = lower_off, ymax = upper_off),
              alpha = 0.35, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(y = est_off, colour = method),           # ← centre = mean ρ̂−ρ
            linewidth = 0.6) +
  scale_fill_manual(values  = c("NI" = "grey70", "INT" = "steelblue"),
                    name = "Estimator") +
  scale_colour_manual(values = c("NI" = "grey40", "INT" = "steelblue"),
                      guide = "none") +
  labs(
    title = title_expr,
    x     = expression(rho),
    y     = expression(mean(CI)-rho)
  ) +
  theme_minimal(base_size = 14)

print(p_band)
#ggsave("fig1_mean_band_vs_rho.pdf", p_band,
#       width = 6.8, height = 4.4)
ggsave("fig1_mean_band_vs_rho_noramlised.pdf", p_band,
       width = 6.8, height = 4.4)
## ------------------------------------------------------------
## 5.4  Figure 2  – CI width & coverage vs n                 ##
## ------------------------------------------------------------
plot2_df <- summ_all[
  eps1 %in% c(0.5, 1.0, 1.5) & eps2 %in% c(0.5, 1.0, 1.5) &
    rho_true == 0.5                                         # fix rho
]

p2_width <- ggplot(plot2_df,
                   aes(x = n, y = ci_len,
                       colour = interaction(eps1, eps2, sep = ","),
                       linetype = method)) +
  geom_line() + geom_point() +
  scale_x_log10() +
  labs(title = bquote("Average CI width vs n  ("*rho*" = 0.5)"),
       x = "n (log-scale)", y = "Average CI length",
       colour = bquote("("*epsilon[1]*","*epsilon[2]*")"), linetype = "Method") +
  theme_minimal(base_size = 14)

p2_cov <- ggplot(plot2_df,
                 aes(x = n, y = coverage,
                     colour = interaction(eps1, eps2, sep = ","),
                     linetype = method)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = 1 - alpha_level, linetype = "dashed") +
  scale_x_log10() +
  labs(title = bquote("Coverage vs n  ("*rho*" = 0.5)"),
       x = "n (log-scale)", y = "Empirical coverage",
       colour = bquote("("*epsilon[1]*","*epsilon[2]*")"), linetype = "Method") +
  theme_minimal(base_size = 14)

plot(p2_cov)
plot(p2_width)
#ggsave("fig2a_ci_width_vs_n.pdf", p2_width, width = 6, height = 4)
#ggsave("fig2b_coverage_vs_n.pdf", p2_cov,  width = 6, height = 4)
ggsave("fig2a_ci_width_vs_n_normalised.pdf", p2_width, width = 6, height = 4)
ggsave("fig2b_coverage_vs_n_normalised.pdf", p2_cov,  width = 6, height = 4)


## ------------------------------------------------------------
## 5.5  Figure 3  – MSE vs n                                 ##
## ------------------------------------------------------------
plot3_df <- summ_all[rho_true == 0.5]

p3 <- ggplot(plot3_df,
             aes(x = n, y = mse,
                 colour = interaction(eps1, eps2, sep = ","),
                 linetype = method)) +
  geom_line() + geom_point() +
  scale_x_log10() + scale_y_log10() +
  labs(title = bquote("MSE of "*hat(rho)*" vs n  ("*rho*" = 0.5)"),
       x = "n (log-scale)", y = "MSE (log-scale)",
       colour =  bquote("("*epsilon[1]*","*epsilon[2]*")"), linetype = "Method") +
  theme_minimal(base_size = 14)

plot(p3)
#ggsave("fig3_mse_vs_n.pdf", p3, width = 6, height = 4)
ggsave("fig3_mse_vs_n_normalised.pdf", p3, width = 6, height = 4)
message("Simulation & figures complete!  PDF files saved to working dir.")

