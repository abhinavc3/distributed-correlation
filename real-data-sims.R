## ============================================================
## HRS (wave 2): BMI vs Age — DP correlation (NI + INT)
## ============================================================

## ---- Libraries ----
library(dplyr)
library(tidyr)
library(ggplot2)
# library(here)  # uncomment if you want reproducible paths via here::here()

## ---- IO ----
# obj <- readRDS(here::here("hrs_long_panel.rds"))
obj <- readRDS("hrs_long_panel.rds")  # keep your path if you're running interactively

## ---- Wave-level missingness summary ----
by_wave <- obj %>%
  mutate(wave = suppressWarnings(as.integer(wave))) %>%
  filter(!is.na(wave)) %>%
  group_by(wave) %>%
  summarise(
    n               = n(),
    missing_age     = sum(is.na(agey_e)),
    missing_bmi     = sum(is.na(bmi)),
    missing_any     = sum(is.na(agey_e) | is.na(bmi)),
    complete_cases  = sum(!is.na(agey_e) & !is.na(bmi)),
    .groups = "drop"
  ) %>%
  arrange(wave) %>%
  mutate(
    pct_missing_age = round(100 * missing_age / n, 1),
    pct_missing_bmi = round(100 * missing_bmi / n, 1),
    pct_missing_any = round(100 * missing_any / n, 1)
  )

# View(by_wave)

## ---- Choose wave 2 and basic plot ----
wave2 <- obj %>%
  filter(wave %in% c(2, "2")) %>%
  transmute(hhidpn, age = agey_e, bmi) %>%
  drop_na(age, bmi)

# Sanity
# nrow(wave2); dplyr::n_distinct(wave2$hhidpn)

# Scatter + smooth
# ggplot(wave2, aes(age, bmi)) +
#   geom_point(alpha = 0.15, size = 0.8) +
#   geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = TRUE) +
#   labs(title = "HRS Wave 2: BMI vs Age", x = "Age (years)", y = "BMI") +
#   theme_minimal()

## ============================================================
## DP building blocks
## ============================================================

## Laplace sampler (central & local DP noise)
rLap <- function(n, scale) {
  u <- runif(n, -0.5, 0.5)
  -scale * sign(u) * log(1 - 2 * abs(u))
}

## DP mean with clipping: sensitivity = (hi - lo)/n
dp_mean <- function(x, lo, hi, eps) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  x_clip <- pmin(pmax(x, lo), hi)
  n <- length(x_clip)
  mean(x_clip) + rLap(1, scale = (hi - lo) / (n * eps))
}

## DP sd via clipped 2nd moment
dp_sd <- function(x, lo, hi, eps1, eps2) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  x_clip <- pmin(pmax(x, lo), hi)
  n <- length(x_clip)
  
  mu_dp <- dp_mean(x_clip, lo, hi, eps1)             # private mean
  m2_dp <- mean(x_clip^2) + rLap(1, (hi^2 - lo^2) / (n * eps2))  # private E[x^2]
  
  sd_dp <- sqrt(max(m2_dp - mu_dp^2, 0))
  list(mean = mu_dp, sd = sd_dp)
}

## Standardize using your private mean/sd and the same clipping bounds
standardize_dp <- function(x, priv, lo, hi, eps = 1e-8) {
  x_clipped <- pmin(pmax(x, lo), hi)
  (x_clipped - priv$mean) / max(priv$sd, eps)
}

standardize_age_bmi <- function(age, bmi,
                                age_priv, bmi_priv,
                                age_lo, age_hi, bmi_lo, bmi_hi,
                                eps = 1e-8) {
  tibble::tibble(
    age_z = standardize_dp(age, age_priv, age_lo, age_hi, eps),
    bmi_z = standardize_dp(bmi, bmi_priv, bmi_lo, bmi_hi, eps)
  )
}

## Symmetric λ for standardized variable from known raw bounds + private mean/sd
lambda_from_priv <- function(lo, hi, priv, eps_sd = 1e-8) {
  sig <- max(priv$sd, eps_sd)
  max(abs((lo - priv$mean) / sig), abs((hi - priv$mean) / sig))
}

## Generic λ(n) if you don’t override (kept for fallback)
lambda_n <- function(n, eta = 1) min(2 * eta * sqrt(log(n)), 2 * sqrt(3))

## ============================================================
## Estimator 1: NI (no interaction) — batched
## ============================================================

correlation_NI_subG <- function(X, Y, eps1, eps2,
                                eta1 = 1, eta2 = 1,
                                alpha = 0.05,
                                lambda_X = NULL, lambda_Y = NULL) {
  ok <- !(is.na(X) | is.na(Y))
  X <- X[ok]; Y <- Y[ok]
  n <- length(X); stopifnot(n == length(Y), n >= 2)
  
  λ1 <- if (!is.null(lambda_X)) lambda_X else lambda_n(n, eta1)
  λ2 <- if (!is.null(lambda_Y)) lambda_Y else lambda_n(n, eta2)
  
  Xc <- pmax(pmin(X,  λ1), -λ1)
  Yc <- pmax(pmin(Y,  λ2), -λ2)
  
  m <- ceiling(8 / (eps1 * eps2)); if (m > n) m <- n
  k <- floor(n / m); if (k < 2) { k <- 2; m <- floor(n / k) }
  idx <- sample.int(n, k * m)  # randomize batches
  
  X_bar <- rowMeans(matrix(Xc[idx], nrow = k, byrow = TRUE))
  Y_bar <- rowMeans(matrix(Yc[idx], nrow = k, byrow = TRUE))
  
  X_tilde <- X_bar + rLap(k, scale = 2 * λ1 / (m * eps1))
  Y_tilde <- Y_bar + rLap(k, scale = 2 * λ2 / (m * eps2))
  
  rho_hat <- (m / k) * sum(X_tilde * Y_tilde)
  
  Tj <- m * X_tilde * Y_tilde
  se <- sd(Tj) / sqrt(k)
  crit <- qnorm(1 - alpha / 2)
  ci <- c(max(rho_hat - crit * se, -1), min(rho_hat + crit * se, 1))
  
  list(rho_hat = rho_hat, ci = ci, k = k, m = m, lambda_X = λ1, lambda_Y = λ2)
}

## ============================================================
## Estimator 2: INT (one-round interaction) + CI via mixture quantile
## ============================================================

# Paper-like defaults for λ if not overridden (fallbacks)
lambda_INT_n <- function(n, eta_s = 1, eta_r = 1, eps_s = 1) {
  lambda_s <- min(2 * eta_s * sqrt(log(n)), 2 * sqrt(3))
  lambda_r <- 5 * max(eta_r, 1) * min(log(n), 6) / (min(eps_s, 1))
  c(lambda_s, lambda_r)  # NOTE: lambda_r not same as paper; used only if no overrides
}

# Quantile of Z + c * (sign)*Exp(1), used for the mixture-quantile CI
mixquant <- function(c, p, nsim = 2000L) {
  x <- rnorm(nsim) + c * rexp(nsim) * (2 * rbinom(nsim, 1, 0.5) - 1)
  sort(x)[ceiling(p * nsim)]
}

# ---- NEW: receiver-λ that accounts for sender's local DP noise ----
# If sender sends S = clip(X, ±λs) + Lap(0, b_s), b_s = 2λs/εs,
# and receiver multiplies by the other variable clipped to ±λo, then
# with prob ≥ 1 - δ per sample: |U| ≤ (λs + b_s log(1/δ)) * λo.
lambda_receiver_from_noise <- function(lambda_sender, lambda_other,
                                       eps_sender, delta_per_sample) {
  b_s <- 2 * lambda_sender / eps_sender
  (lambda_sender + b_s * log(1 / delta_per_sample)) * lambda_other
}

ci_INT_subG <- function(X, Y, eps1, eps2,
                        eta1 = 1, eta2 = 1,
                        alpha = 0.05,
                        mode  = c("auto","normal","laplace"),
                        # OPTIONAL: bounds (standardized variables recommended)
                        lambda_sender   = NULL,   # bound for sender variable (±)
                        lambda_other    = NULL,   # bound for non-sender variable (±)
                        lambda_receiver = NULL,   # bound for product; computed if NULL
                        delta_clip      = NULL    # per-sample tail prob; default 1/n
) {
  # 0) clean & pairwise complete
  ok <- !(is.na(X) | is.na(Y))
  X <- X[ok]; Y <- Y[ok]
  n <- length(X); stopifnot(n == length(Y), n >= 2)
  
  # 1) roles: larger epsilon is sender
  sender_is_X <- (eps1 >= eps2)
  eps_s <- if (sender_is_X) eps1 else eps2
  eps_r <- if (sender_is_X) eps2 else eps1
  eta_s <- if (sender_is_X) eta1 else eta2
  eta_r <- if (sender_is_X) eta2 else eta1
  
  # 2) default per-sample tail level
  if (is.null(delta_clip)) delta_clip <- 1 / n  # or, e.g., 0.01/n
  
  # 3) get lambdas (use fallbacks only if overrides not provided)
  if (is.null(lambda_sender) || is.null(lambda_other)) {
    λ <- lambda_INT_n(n, eta_s = eta_s, eta_r = eta_r, eps_s = eps_s)
    lambda_sender <- if (is.null(lambda_sender)) λ[1] else lambda_sender
    if (is.null(lambda_other)) {
      lambda_other <- lambda_n(n, if (sender_is_X) eta2 else eta1)
    }
  }
  
  # 4) receiver bound that accounts for sender noise (if not provided)
  if (is.null(lambda_receiver)) {
    lambda_receiver <- lambda_receiver_from_noise(
      lambda_sender    = lambda_sender,
      lambda_other     = lambda_other,
      eps_sender       = eps_s,
      delta_per_sample = delta_clip
    )
  }
  
  # 5) sender clips & adds local DP noise; receiver clips other var and multiplies
  if (sender_is_X) {
    Xc <- pmax(pmin(X,  lambda_sender), -lambda_sender)
    Yb <- pmax(pmin(Y,  lambda_other ), -lambda_other )
    U  <- (Xc + rLap(n, scale = 2 * lambda_sender / eps_s)) * Yb
  } else {
    Yc <- pmax(pmin(Y,  lambda_sender), -lambda_sender)
    Xb <- pmax(pmin(X,  lambda_other ), -lambda_other )
    U  <- (Yc + rLap(n, scale = 2 * lambda_sender / eps_s)) * Xb
  }
  
  # 6) receiver clips product and adds one central-DP Laplace draw to the MEAN
  Uc <- pmax(pmin(U, lambda_receiver), -lambda_receiver)
  rho_hat <- mean(Uc) + rLap(1, scale = 2 * lambda_receiver / (n * eps_r))
  
  # 7) CI via mixture-quantile (Z + c*Lap), using sampling SE only
  sdUc <- sd(Uc)
  if (sdUc == 0) {
    width <- qnorm(1 - alpha/2) * sqrt(2) * (2 * lambda_receiver / (n * eps_r))
  } else {
    cstar <- (2 * lambda_receiver) / (sqrt(n) * sdUc * eps_r)
    width <- mixquant(cstar, 1 - alpha/2) * (sdUc / sqrt(n))
  }
  ci <- c(max(rho_hat - width, -1), min(rho_hat + width, 1))
  
  list(rho_hat = rho_hat,
       ci      = ci,
       roles   = if (sender_is_X) "X→Y" else "Y→X",
       lambda_sender   = lambda_sender,
       lambda_other    = lambda_other,
       lambda_receiver = lambda_receiver,
       delta_clip      = delta_clip)
}


## ============================================================
## Run: private standardization + λ overrides from known ranges
## ============================================================

# Clipping bounds in raw units
age_lo <- 45; age_hi <- 90
bmi_lo <- 15; bmi_hi <- 35

# Privacy budgets
eps_mean  <- 0.10    # mean
eps_m2    <- 0.10    # second moment
eps_corr <- 2.00
eps_corr1 <- eps_corr    # for NI sender X
eps_corr2 <- eps_corr    # for NI sender Y (same if symmetric)
eps_int_X <- eps_corr    # INT: epsilon for X side
eps_int_Y <- eps_corr    # INT: epsilon for Y side

# Private mean/sd (central DP)
age_priv <- dp_sd(wave2$age, age_lo, age_hi, eps_mean, eps_m2)
bmi_priv <- dp_sd(wave2$bmi, bmi_lo, bmi_hi, eps_mean, eps_m2)

# Standardize + drop NAs (should be none after our drop above)
wave2_std <- standardize_age_bmi(
  age = wave2$age, bmi = wave2$bmi,
  age_priv = age_priv, bmi_priv = bmi_priv,
  age_lo = age_lo, age_hi = age_hi,
  bmi_lo = bmi_lo, bmi_hi = bmi_hi
)
wave2_std_clean <- drop_na(wave2_std, age_z, bmi_z)

# λ from private standardization + known raw ranges
lambda_age_z <- lambda_from_priv(age_lo, age_hi, age_priv)
lambda_bmi_z <- lambda_from_priv(bmi_lo, bmi_hi, bmi_priv)

# ---- NI (batched) on standardized data with λ overrides ----
set.seed(231)
res_NI <- correlation_NI_subG(
  X = wave2_std_clean$age_z,
  Y = wave2_std_clean$bmi_z,
  eps1 = eps_corr1, eps2 = eps_corr2,
  lambda_X = lambda_age_z, lambda_Y = lambda_bmi_z
)

# ---- INT (one-round) example: let AGE be the sender ----
# Per standardized-variable bounds:
lambda_age_z <- lambda_from_priv(age_lo, age_hi, age_priv)  # sender bound if AGE sends
lambda_bmi_z <- lambda_from_priv(bmi_lo, bmi_hi, bmi_priv)  # other side's bound

# Receiver λ that accounts for AGE's local DP noise.
# Use per-sample δ; union-bound choice δ = 1/n is a good start.
delta_clip <- 1 / nrow(wave2_std_clean)
lambda_receiver_ageSends <- lambda_receiver_from_noise(
  lambda_sender    = lambda_age_z,
  lambda_other     = lambda_bmi_z,
  eps_sender       = eps_int_X,        # sender's ε
  delta_per_sample = delta_clip
)

set.seed(322)
res_INT <- ci_INT_subG(
  X = wave2_std_clean$age_z,
  Y = wave2_std_clean$bmi_z,
  eps1 = eps_int_X,  # AGE epsilon (sender)
  eps2 = eps_int_Y,  # BMI epsilon (receiver)
  lambda_sender   = lambda_age_z,
  lambda_other    = lambda_bmi_z,
  lambda_receiver = lambda_receiver_ageSends,
  delta_clip      = delta_clip
)

## ---- Results ----
cat("\n--- Private standardization ---\n")
print(age_priv); print(bmi_priv)

cat("\n--- NI estimator (standardized, clipped) ---\n")
print(res_NI)

cat("\n--- INT estimator (AGE→BMI) ---\n")
print(res_INT)




library(dplyr)
library(tidyr)
library(ggplot2)

## ---------------------------------------
## Settings
## ---------------------------------------
eps_grid <- seq(0.25, 2.5, by = 0.1)  # x-axis
R <- 200                              # replications per epsilon (tune as you like)

# Non-private baseline (on standardized, clipped data)
rho_np <- cor(wave2_std_clean$age_z, wave2_std_clean$bmi_z, use = "complete.obs")

## ---------------------------------------
## Single-run helpers (use your functions)
## ---------------------------------------

run_NI_once <- function(eps_corr, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  res <- correlation_NI_subG(
    X = wave2_std_clean$age_z,
    Y = wave2_std_clean$bmi_z,
    eps1 = eps_corr,
    eps2 = eps_corr,
    lambda_X = lambda_age_z,
    lambda_Y = lambda_bmi_z
  )
  tibble(
    method   = "NI",
    eps_corr = eps_corr,
    rho_hat  = res$rho_hat,
    ci_low   = res$ci[1],
    ci_high  = res$ci[2]
  )
}

run_INT_once <- function(eps_corr, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(wave2_std_clean)
  delta_clip <- 1 / n  # union-bound-ish choice
  
  # Receiver λ that accounts for sender's local DP noise (AGE is sender here)
  lambda_receiver <- lambda_receiver_from_noise(
    lambda_sender    = lambda_age_z,
    lambda_other     = lambda_bmi_z,
    eps_sender       = eps_corr,
    delta_per_sample = delta_clip
  )
  
  res <- ci_INT_subG(
    X = wave2_std_clean$age_z,  # AGE (sender)
    Y = wave2_std_clean$bmi_z,
    eps1 = eps_corr,  # sender ε
    eps2 = eps_corr,  # receiver ε (same here; change if you want asymmetry)
    lambda_sender   = lambda_age_z,
    lambda_other    = lambda_bmi_z,
    lambda_receiver = lambda_receiver,
    delta_clip      = delta_clip
  )
  tibble(
    method   = "INT",
    eps_corr = eps_corr,
    rho_hat  = res$rho_hat,
    ci_low   = res$ci[1],
    ci_high  = res$ci[2]
  )
}

## ---------------------------------------
## Run many reps and summarize by epsilon
## ---------------------------------------

# Build a grid of (eps, rep)
grid <- tidyr::expand_grid(eps_corr = eps_grid, rep = seq_len(R))

# NI runs
ni_runs <- grid |>
  rowwise() |>
  do(run_NI_once(.$eps_corr, seed = 10 + 37*.$rep + 1000*which(eps_grid == .$eps_corr))) |>
  ungroup()

ni_mean <- ni_runs |>
  group_by(eps_corr) |>
  summarise(
    method      = "NI",
    rho_hat_mean = mean(rho_hat),
    ci_low_mean  = mean(ci_low),
    ci_high_mean = mean(ci_high),
    # Optional spread (for ribbons / error bars if you want)
    ci_low_q10   = quantile(ci_low, 0.10),
    ci_high_q90  = quantile(ci_high, 0.90),
    .groups = "drop"
  )

# INT runs
int_runs <- grid |>
  rowwise() |>
  do(run_INT_once(.$eps_corr, seed = 20 + 41*.$rep + 1000*which(eps_grid == .$eps_corr))) |>
  ungroup()

int_mean <- int_runs |>
  group_by(eps_corr) |>
  summarise(
    method       = "INT",
    rho_hat_mean = mean(rho_hat),
    ci_low_mean  = mean(ci_low),
    ci_high_mean = mean(ci_high),
    ci_low_q10   = quantile(ci_low, 0.10),
    ci_high_q90  = quantile(ci_high, 0.90),
    .groups = "drop"
  )

## ---------------------------------------
## Plots: mean CI per epsilon + mean rho_hat
## ---------------------------------------

library(dplyr)
library(ggplot2)
library(patchwork)   # install.packages("patchwork") if needed

# 1) Compute mean(CI) per epsilon for both methods
ni_mean <- ni_mean %>% mutate(mean_ci = (ci_low_mean + ci_high_mean) / 2)
int_mean <- int_mean %>% mutate(mean_ci = (ci_low_mean + ci_high_mean) / 2)

# 2) Shared y-axis limits (include CIs, non-private rho, and zero)
y_all <- c(ni_mean$ci_low_mean, ni_mean$ci_high_mean,
           int_mean$ci_low_mean, int_mean$ci_high_mean,
           rho_np, 0)
y_pad <- 0.02 * diff(range(y_all, na.rm = TRUE))
y_lim <- range(y_all, na.rm = TRUE) + c(-y_pad, y_pad)

# 3) A consistent, clean theme
paper_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6))
  )

# 4) Panels
p_ni_mean <- ggplot(ni_mean, aes(eps_corr)) +
  geom_errorbar(aes(ymin = ci_low_mean, ymax = ci_high_mean), width = 0.08) +
  geom_point(aes(y = mean_ci), size = 1.8) +
  geom_hline(yintercept = rho_np, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "red") +
  coord_cartesian(ylim = y_lim) +
  labs(
    title = "Non-interactive",
    x = expression(epsilon[corr]),
    y = expression("mean(CI) for " * rho)
  ) +
  paper_theme

p_int_mean <- ggplot(int_mean, aes(eps_corr)) +
  geom_errorbar(aes(ymin = ci_low_mean, ymax = ci_high_mean), width = 0.08) +
  geom_point(aes(y = mean_ci), size = 1.8) +
  geom_hline(yintercept = rho_np, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "red") +
  coord_cartesian(ylim = y_lim) +
  labs(
    title = "Interactive",
    x = expression(epsilon[corr]),
    y = expression("mean(CI) for " * rho)
  ) +
  paper_theme

# 5) Side-by-side with shared y-axis
p_combined <- p_ni_mean + p_int_mean + plot_layout(ncol = 2)
p_combined

