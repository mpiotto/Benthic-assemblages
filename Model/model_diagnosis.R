#=============================================================================
# Posterior Predictive Checks for the JSDM
#-----------------------------------------------------------------------------
# 1. Convergence diagnostics (Rhat and n_eff) on parameters with priors.
# 2. Posterior predictive checks for occupancy and conditional cover.
#=============================================================================

# ---- Setup ----------------------------------------------------------------

rm(list = ls())

library(rstan)
library(ggplot2)
library(dplyr)

# Reusable plot theme
theme_clean <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      axis.text        = element_text(colour = "black", size = 12),
      axis.line        = element_line(colour = "black"),
      axis.title       = element_text(size = 16),
      legend.title     = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Load fitted model and input data
load("inputdata.RData")
load("model_fit.RData")

fit_summary <- summary(fit, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))$summary
draws       <- rstan::extract(fit)


# ---- 1. Convergence: Rhat and n_eff ---------------------------------------
# Diagnostics are only meaningful on parameters with priors (not derived ones).

priors_named <- c("tau", "rho", "phi", "rhosq", "etasq", "delta")
idx_priors   <- which(
  rownames(fit_summary) %in% priors_named |
    grepl("^z(\\[|$)", rownames(fit_summary))
)

diag_df <- data.frame(
  Rhat = fit_summary[idx_priors, "Rhat"],
  neff = fit_summary[idx_priors, "n_eff"]
)

cat(sprintf("n_eff  : min = %.0f   max = %.0f\n",
            min(diag_df$neff, na.rm = TRUE),
            max(diag_df$neff, na.rm = TRUE)))
cat(sprintf("Rhat   : min = %.3f  max = %.3f\n",
            min(diag_df$Rhat, na.rm = TRUE),
            max(diag_df$Rhat, na.rm = TRUE)))

# Fig. SA3
p_rhat <- ggplot(diag_df, aes(x = Rhat)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.65, size = 0.4) +
  labs(x = "Rhat", y = "Density") +
  theme_clean()

p_neff <- ggplot(diag_df, aes(x = neff)) +
  geom_histogram(binwidth = 4000, colour = "black",
                 fill = "grey", alpha = 0.65, size = 0.4) +
  labs(x = "Mean effective sample size", y = "Frequency") +
  theme_clean()

print(p_rhat)
print(p_neff)


# ---- 2. Posterior predictive simulation -----------------------------------
# Coefficient layout (16 per species). 1-8: presence (psi). 9-16: cover (mu).
# 
#  pres   cover        parameter
#   1 and  9 : intercept (~15 m, 1994)
#   2 and 10 : intercept differential, 1998
#   3 and 11 : intercept differential, 2010
#   4 and 12 : intercept differential, 2022
#   5 and 13 : slope w.r.t. bathymetry
#   6 and 14 : slope differential, 1998
#   7 and 15 : slope differential, 2010
#   8 and 16 : slope differential, 2022

K_pres  <- 8
K_cov   <- 8
K_total <- K_pres + K_cov

# Reshape draws$betas, which is [n_draws, nsp * K_total] in column-major
# order (Stan stacks columns of b: all species for coef 1, then coef 2, ...),
# into a tidy 3-D array [draw, species, coef].
n_draws   <- nrow(draws$betas)
betas_arr <- array(draws$betas, dim = c(n_draws, nsp, K_total))

# Subsample posterior draws to keep simulation manageable.
# (We could use all draws; sampling just thins the posterior cheaply.)
n_iter <- 1000
iter   <- sample.int(n_draws, n_iter, replace = FALSE)

betas_iter  <- betas_arr[iter, , ]                         # [n_iter, nsp, K]
phi_iter    <- draws$phi[iter]                             # [n_iter]
delta_tr    <- draws$delta_tr[iter, ]                      # [n_iter, n_T]
delta_tr_sp <- aperm(draws$delta_tr_sp[iter, , ],          # -> [n_T, nsp, n_iter]
                     c(2, 3, 1))

# Design matrices for presence and cover (currently identical).
X_psi <- as.matrix(X); rownames(X_psi) <- NULL
X_mu  <- X_psi

# Storage: predicted occupancy and predicted conditional mean cover.
occ_pred   <- matrix(NA_real_, n_iter, nsp)
cover_pred <- matrix(NA_real_, n_iter, nsp)

for (ii in seq_len(n_iter)) {
  for (jj in seq_len(nsp)) {
    
    # --- Presence component ---
    beta_pres <- betas_iter[ii, jj, 1:K_pres]
    ranef     <- delta_tr[ii, transect_id] +
      delta_tr_sp[transect_id, jj, ii]
    psi       <- plogis(X_psi %*% beta_pres + ranef)
    z         <- rbinom(length(psi), size = 1, prob = psi)
    occ_pred[ii, jj] <- mean(z)
    
    # --- Cover component (only meaningful where z == 1) ---
    beta_cov <- betas_iter[ii, jj, (K_pres + 1):K_total]
    mu       <- plogis(X_mu %*% beta_cov)
    p_       <- mu * phi_iter[ii]
    q_       <- phi_iter[ii] - p_
    y_sim    <- rbeta(length(mu), p_, q_)
    cover_pred[ii, jj] <- if (any(z == 1)) mean(y_sim[z == 1]) else NA_real_
  }
}


# ---- 3. Predicted vs observed plots ---------------------------------------

# Observed quantities (uncomment to include observed points in plots)
# obs_occ <- apply(Y, 2, function(x) mean(x > 0))
# obs_cov <- apply(Y, 2, function(x) mean(x[x > 0]))

summarise_pred <- function(M) {
  data.frame(
    pred = colMeans(M, na.rm = TRUE),
    L    = apply(M, 2, quantile, probs = 0.025, na.rm = TRUE),
    U    = apply(M, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

df_occ <- cbind(species = species, summarise_pred(occ_pred))   # , obs = obs_occ
df_cov <- cbind(species = species, summarise_pred(cover_pred)) # , obs = obs_cov

# Fig. SA4: predicted in green.
# Observed points in purple -- uncomment geom_point and obs columns above to show.
ppc_plot <- function(df, ylab) {
  ggplot(df, aes(x = species)) +
    geom_errorbar(aes(ymin = L, ymax = U),
                  width = 0.05, size = 1, colour = "#7AD151FF") +
    geom_point(aes(y = pred), size = 4, colour = "#7AD151FF") +
    # geom_point(aes(y = obs), size = 4, colour = "#414487FF", alpha = 0.8) +
    labs(x = "Species", y = ylab) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90,
                                     margin = margin(t = 1),
                                     size = 10),
          panel.grid.minor = element_blank())
}

print(ppc_plot(df_occ, "Occupancy"))
print(ppc_plot(df_cov, "Conditional cover"))
        

