lambda_n <- function(n, eta = 1)  min(2 * eta * sqrt(log(n)),2*sqrt(3))   # same as paper

lambda_INT_n <- function(n, eta_s = 1, eta_r = 1, eps_s = 1){
  lambda_s =  min(2 * eta_s * sqrt(log(n)),2*sqrt(3))
  lambda_r = 5 * max(eta_r,1) * min(log(n),6)/ (min(eps_s, 1))
  return (c(lambda_s, lambda_r))
}    # lambda_r NOT same as paper
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

############################################################################
##  NI  clipped-batch estimator  and CI (sub-Gaussian)                            ##
############################################################################
correlation_NI_subG  <- function(X, Y, eps1, eps2,
                                 eta1 = 1, eta2 = 1,
                                 alpha = 0.05) {
  n <- length(X);  stopifnot(n == length(Y))
  
  λ1 <- lambda_n(n, eta1)
  λ2 <- lambda_n(n, eta2)
  
  Xc <- pmax(pmin(X, λ1), -λ1)      # 1) clip once
  Yc <- pmax(pmin(Y, λ2), -λ2)
  
  ## ---------- batch design (same m, k) -------------------------------
  m <- ceiling(8 / (eps1 * eps2));  if (m > n) m <- n
  k <- floor(n / m);  stopifnot(k >= 1)
  
  idx <- seq_len(k * m)
  X_mat <- matrix(Xc[idx], nrow = k, byrow = TRUE)
  Y_mat <- matrix(Yc[idx], nrow = k, byrow = TRUE)
  
  X_bar <- rowMeans(X_mat)
  Y_bar <- rowMeans(Y_mat)
  
  ## add Laplace noise – sens = 2λ/m
  X_tilde <- X_bar + rLap(k, scale = 2 * λ1 / (m * eps1))
  Y_tilde <- Y_bar + rLap(k, scale = 2 * λ2 / (m * eps2))
  
  eta_hat <- (m / k) * sum(X_tilde * Y_tilde)   # unbiased for ρ
  rho_hat <- eta_hat                            # (no sine link)
  
  ## ---------- CI ------------------------------------------------------
  Tj   <- m * X_tilde * Y_tilde                 # k batch products
  se   <- sd(Tj) / sqrt(k)                     # asymptotic s.e.
  crit <- qnorm(1 - alpha/2)
  ci   <- c(max(rho_hat - crit*se, -1),
            min(rho_hat + crit*se,  1))
  
  list(rho_hat = rho_hat, ci = ci)
}

############################################################################
##  INT clipped estimator  and CI (sub-Gaussian)                          ##
############################################################################
ci_INT_subG <- function(X, Y, eps1, eps2,
                        eta1 = 1, eta2 = 1,
                        alpha = 0.05,
                        mode  = c("auto","normal","laplace")) {
  
  
  n <- length(X);  stopifnot(n == length(Y))
  
  ## choose roles by larger ε
  sender_is_X <- (eps1 >= eps2)
  eps_s <- if (sender_is_X) eps1 else eps2
  eps_r <- if (sender_is_X) eps2 else eps1
  
  eta_s <- if (sender_is_X) eta1 else eta2
  eta_r <- if (sender_is_X) eta2 else eta1
  
  λ <- lambda_INT_n(n, eta_s = eta_s, eta_r = eta_r, eps_s = eps_s)
  λs <- λ[1]
  λr <- λ[2]
  
  if(sender_is_X){
    Xc <- pmax(pmin(X, λs), -λs)
    U <- (Xc + rLap(n, scale = 2 * λs / (eps_s))) * Y
    Uc <- pmax(pmin(U, λr), -λr)
    rho_hat = mean(Uc) + rLap(1, scale = 2 * λr / (n * eps_r))
  } else{
    Yc <- pmax(pmin(Y, λs), -λs)
    U <- (Yc + rLap(n, scale = 2 * λs / (eps_s))) * X
    Uc <- pmax(pmin(U, λr), -λr)
    rho_hat = mean(Uc) + rLap(1, scale = 2 * λr / (n * eps_r))
  }
  
  se_norm = sqrt(sd(Uc)^2 + 2 *(2 * λr / (n * eps_r))^2)
  cstar <- 2/(sqrt(n)*sd(Uc)*eps_r)
  width <- mixquant(cstar, 1-alpha/2) * se_norm/sqrt(n)
  ci   <- c(max(rho_hat - width, -1),
            min(rho_hat + width,  1))
  list(rho_hat = rho_hat,
       ci      = ci,
       mode    = mode,
       roles   = if (sender_is_X) "X→Y" else "Y→X")
}



############################################################
##  DGP: Correlated Gaussian mixture                      ##
############################################################
gen_mix_gaussian <- function(n, rho,
                             mu0    = c(0, 0), sigma0 = c(1, 1),   # comp-0
                             mu1    = c(3, 3), sigma1 = c(2, 0.5), # comp-1
                             pi_mix = 0.5) {
  
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
  
  Σ0 <- matrix(c(sigma0[1]^2, sigma0[1]*sigma0[2]*rho,
                 sigma0[1]*sigma0[2]*rho, sigma0[2]^2), 2)
  Σ1 <- matrix(c(sigma1[1]^2, sigma1[1]*sigma1[2]*rho,
                 sigma1[1]*sigma1[2]*rho, sigma1[2]^2), 2)
  
  labels <- rbinom(n, 1, pi_mix)
  n0 <- sum(labels == 0);   n1 <- n - n0
  
  out <- rbind(
    MASS::mvrnorm(n0, mu = mu0, Sigma = Σ0),
    MASS::mvrnorm(n1, mu = mu1, Sigma = Σ1)
  )
  out = out[sample.int(n), ]# shuffle rows
  pmax(pmin(out, 1), -1)
}

################################################################
##  DGP: bounded common-factor (mean 0, var 1, corr = ρ)      ##
################################################################
gen_bounded_factor <- function(n, rho) {
 # stopifnot(abs(rho) <= 1, L > 0)
  #
  #  X = U + E1 ,  Y = U + E2
  #  Var(U)  = ρ     ( ⇒  U ~ Unif[-cU,cU]  with  cU = sqrt(3ρ) )
  #  Var(Ei) = 1-ρ   ( ⇒  Ei ~ Unif[-cE,cE] with  cE = sqrt(3(1-ρ)) )
  #
  cU <- sqrt(3 * rho)
  cE <- sqrt(3 * (1 - rho))
  U  <- runif(n, -cU,  cU)
  E1 <- runif(n, -cE,  cE)
  E2 <- runif(n, -cE,  cE)
  cbind(U + E1, U + E2)
}

################################################################
##  Simulator  (works for NI/INT Gaussian-sign  *or* sub-G)   ##
################################################################
run_sim_one <- function(n, rho,
                        eps1, eps2,
                        dgp_fun  = gen_bounded_factor,  # ← default DGP
                        dgp_args = list(),              # extra args for DGP
                        B        = 1000,
                        alpha    = 0.05,
                        use_subG = TRUE,                # TRUE = sub-G estimators
                        ci_mode  = "auto",
                        seed     = 2025L) {
  
  set.seed(seed)
  out <- data.frame(repl = seq_len(B),
                    ni_hat = NA, ni_low = NA, ni_up = NA,
                    int_hat = NA, int_low = NA, int_up = NA)
  
  for (b in seq_len(B)) {
    XY <- do.call(dgp_fun, c(list(n = n, rho = rho), dgp_args))
    X  <- XY[, 1];  Y <- XY[, 2]
    
    ## ---------- non-interactive ---------------------------
    if (use_subG) {
      ni <- correlation_NI_subG(X, Y, eps1, eps2, alpha = alpha)
    } else {
      ni <- ci_NI_signbatch(X, Y, eps1, eps2,
                            alpha = alpha, normalise = TRUE)
    }
    out$ni_hat[b] <- ni$rho_hat
    out$ni_low[b] <- ni$ci[1];  out$ni_up[b] <- ni$ci[2]
    
    ## ---------- interactive -------------------------------
    if (use_subG) {
      int <- ci_INT_subG(X, Y, eps1, eps2, alpha = alpha)
    } else {
      int <- ci_INT_signflip(X, Y, eps1, eps2,
                             alpha = alpha, mode = ci_mode,
                             normalise = TRUE)
    }
    out$int_hat[b] <- int$rho_hat
    out$int_low[b] <- int$ci[1]; out$int_up[b] <- int$ci[2]
  }
  
  ## convenience metrics -----------------------------------
  out$ni_se2  <- (out$ni_hat  - rho)^2
  out$int_se2 <- (out$int_hat - rho)^2
  out$ni_cover  <- rho >= out$ni_low  & rho <= out$ni_up
  out$int_cover <- rho >= out$int_low & rho <= out$int_up
  out$ni_ci_len <- out$ni_up  - out$ni_low
  out$int_ci_len<- out$int_up - out$int_low
  
  summarise <- function(est,se2,cov,len)
    c(mse = mean(se2), bias = mean(est) - rho, var = var(est),
      coverage = mean(cov), ci_length = mean(len))
  
  summary <- rbind(
    NI  = summarise(out$ni_hat,  out$ni_se2,
                    out$ni_cover,  out$ni_ci_len),
    INT = summarise(out$int_hat, out$int_se2,
                    out$int_cover, out$int_ci_len)
  )
  summary <- data.frame(method = rownames(summary), summary,
                        row.names = NULL)
  
  list(detail = out, summary = summary)
}

res <- run_sim_one(
  n   = 5500,
  rho = 0.6,
  eps1 = 5, eps2 = 1,
  B = 500,   
  use_subG  = TRUE          # FALSE → Gaussian-sign pipeline
)

print(res$summary)
head(res$detail)



################################################################
##  5.  GRID-WISE SIMULATION  – sub-Gaussian mixture DGP      ##
##      – NI clipped-batch  &  INT clipped sign-flip          ##
################################################################

## ------------------------------------------------------------
## 5.1  Design grid & DGP parameters                          ##
## ------------------------------------------------------------
n_grid   <- c(2500, 4000, 6000, 9000, 12000)
rho_grid <- c(0, 0.15, 0.3, 0.4, 0.5, 0.65, 0.8, 0.9)

eps_pairs <- list(
  c(0.5, 0.5),
  c(1.0, 1.0),
  c(1.5, 0.5)      # asymmetric
)



B_per_setting <-  250     # ← raise for camera-ready
alpha_level   <- 0.05
ci_mode_choice<- "auto"

## ------------------------------------------------------------
## 5.2  Run the grid (parallel where possible)                ##
## ------------------------------------------------------------
expand_grid <- function(...) expand.grid(..., KEEP.OUT.ATTRS = FALSE)

design_df <- expand_grid(
  n       = n_grid,
  rho     = rho_grid,
  eps_idx = seq_along(eps_pairs)
)

n_cores <- max(1L, parallel::detectCores() - 1L)

run_row <- function(i) {
  row  <- design_df[i, ]
  eps1 <- eps_pairs[[row$eps_idx]][1]
  eps2 <- eps_pairs[[row$eps_idx]][2]
  
  run_sim_one(
    n       = row$n,
    rho     = row$rho,
    eps1    = eps1,
    eps2    = eps2,
    dgp_fun = gen_bounded_factor,
    use_subG= TRUE,               # << use sub-Gaussian estimators
    B       = B_per_setting,
    alpha   = alpha_level,
    seed    = 1e6 + i
  )
}

if (.Platform$OS.type == "windows") {
  grid_results <- lapply(seq_len(nrow(design_df)), run_row)
} else {
  grid_results <- parallel::mclapply(
    seq_len(nrow(design_df)), run_row, mc.cores = n_cores)
}

## ------------------------------------------------------------
## 5.2.1  Merge replicate-level data                          ##
## ------------------------------------------------------------
library(data.table)

for (i in seq_len(nrow(design_df))) {
  meta <- design_df[i, ]
  det  <- grid_results[[i]]$detail
  det$n        <- meta$n
  det$rho_true <- meta$rho
  det$eps1     <- eps_pairs[[meta$eps_idx]][1]
  det$eps2     <- eps_pairs[[meta$eps_idx]][2]
  grid_results[[i]]$detail <- det
}

detail_all <- rbindlist(lapply(grid_results, `[[`, "detail"))
saveRDS(detail_all, "subG_sim_detail_all.rds")

## ------------------------------------------------------------
## 5.2.2  Summary tables (NI vs INT)                          ##
## ------------------------------------------------------------
summ_INT <- detail_all[
  , .(mse = mean(int_se2),
      bias = mean(int_hat) - mean(rho_true),
      coverage = mean(int_cover),
      ci_len   = mean(int_ci_len)),
  by = .(n, rho_true, eps1, eps2)
][, method := "INT"]

summ_NI <- detail_all[
  , .(mse = mean(ni_se2),
      bias = mean(ni_hat) - mean(rho_true),
      coverage = mean(ni_cover),
      ci_len   = mean(ni_ci_len)),
  by = .(n, rho_true, eps1, eps2)
][, method := "NI"]

summ_all <- rbindlist(list(summ_NI, summ_INT))

## ------------------------------------------------------------
## 5.3  Figure 1 – mean offset bands vs ρ                     ##
## ------------------------------------------------------------
library(ggplot2)

chosen_n   <- 6000
chosen_eps <- c(1.5, 0.5)

slice <- detail_all[
  n == chosen_n & eps1 == chosen_eps[1] & eps2 == chosen_eps[2]
]

band_dt <- rbind(
  slice[, .(rho = rho_true, lower_off = mean(ni_low  - rho_true),
            upper_off = mean(ni_up   - rho_true),
            est_off   = mean(ni_hat - rho_true),
            method = "NI"),  by = rho_true],
  slice[, .(rho = rho_true, lower_off = mean(int_low - rho_true),
            upper_off = mean(int_up  - rho_true),
            est_off   = mean(int_hat - rho_true),
            method = "INT"), by = rho_true])

band_dt[, method := factor(method, levels = c("NI","INT"))]

title_expr <- bquote("Mean CI offset bands — n = "*.(chosen_n)*
                       ",  "*epsilon[1]*" = "*.(chosen_eps[1])*
                       ",  "*epsilon[2]*" = "*.(chosen_eps[2]))

p_band <- ggplot(band_dt, aes(x = rho, fill = method)) +
  geom_ribbon(aes(ymin = lower_off, ymax = upper_off),
              alpha = 0.35, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(y = est_off, colour = method), linewidth = 0.65) +
  scale_fill_manual(values  = c("NI" = "grey70", "INT" = "steelblue"),
                    name = "Estimator") +
  scale_colour_manual(values = c("NI" = "grey35", "INT" = "steelblue"),
                      guide = "none") +
  labs(title = title_expr,
       x = expression(rho),
       y = expression(mean(CI)-rho)) +
  theme_minimal(base_size = 14)

plot(p_band)
ggsave("subG_fig1_mean_band.pdf", p_band, width = 6.8, height = 4.4)

## ------------------------------------------------------------
## 5.4  Figure 2 – CI width & coverage vs n                  ##
## ------------------------------------------------------------
plot2_df <- summ_all[
  eps1 %in% c(0.5,1,1.5) & eps2 %in% c(0.5,1,1.5) & rho_true == 0.5]

p2_width <- ggplot(plot2_df,
                   aes(x = n, y = ci_len,
                       colour = interaction(eps1, eps2, sep = ","),
                       linetype = method)) +
  geom_line() + geom_point() + scale_x_log10() +
  labs(title = bquote("Average CI width vs n  ("*rho*" = 0.5)"),
       x = "n (log-scale)", y = "Average CI length",
       colour = bquote("("*epsilon[1]*","*epsilon[2]*")"),
       linetype = "Method") +
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
       colour = bquote("("*epsilon[1]*","*epsilon[2]*")"),
       linetype = "Method") +
  theme_minimal(base_size = 14)
plot(p2_cov)
plot(p2_width)
ggsave("subG_fig2a_width.pdf", p2_width, width = 6, height = 4)
ggsave("subG_fig2b_cov.pdf",   p2_cov,   width = 6, height = 4)

## ------------------------------------------------------------
## 5.5  Figure 3 – MSE vs n                                  ##
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
       colour = bquote("("*epsilon[1]*","*epsilon[2]*")"),
       linetype = "Method") +
  theme_minimal(base_size = 14)

plot(p3)
ggsave("subG_fig3_mse.pdf", p3, width = 6, height = 4)

message("sub-Gaussian simulation complete — figures saved.")
