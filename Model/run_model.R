library(tidyverse)
library(rstan)

# Design stan model -------------------------------------------------------

rm(list= ls())

# Data --------------------------------------------------------------------

input_fn <- "inputdata.RData"
load(input_fn)

# Prepare data ------------------------------------------------------------

# List with vectors indicating which rows have presence by genus
pres_rows <- lapply(1:ncol(Y), function(i) {
  unname(which(Y[, i] > 0))
})

# A small fraction of the data has cover > 0. To reduce computation time, the
# mean cover will be computed only at the rows of each genus where cover is
# > 0.
# As Stan does not allow for rugged arrays, we need a vector indicating which
# rows to evaluate for each genus. This will be done using two vectors that
# indicate where each genus starts and ends.

pres_n <- lapply(pres_rows, length) |> unlist()
end <- cumsum(pres_n)
begin <- end - pres_n + 1
rows_cover <- do.call("c", pres_rows) # vectorized list

# Transect design matrix
X_transect <- matrix(0, nrow(Y), max(transect_id))
for(i in 1:max(transect_id)) {
  X_transect[transect_id == i, i] <- 1
}
# View(cbind(X_transect, transect_id)) # OK

# binary presence matrix
Y_pres <- Y
Y_pres[Y > 0] <- 1
# unique(as.vector(Y_pres))
# unique(as.vector(Y)) ## OK

# Stan data ---------------------------------------------------------------

stan_dat <- list(
  N = nrow(Y),
  Nt = max(transect_id),
  X_transect = X_transect,
  Kpres = ncol(X),
  Kcov = ncol(X),
  Xpres = X,
  Xcov = X,
  K = ncol(X)*2,
  N_J = ngen,
  L_J = ncol(TT),
  TT = TT,
  C = phyl,
  Y_pres = Y_pres,
  Y_cov = Y,
  prior_sd_z = 2,

  # data to evaluate cover likelihood
  Npres = sum(Y > 0),
  rows_cover = rows_cover,
  begin = begin,
  pres_length = pres_n
)

# Stan settings -----------------------------------------------------------

pars <- c("z", "betas", "phi",
          "Omega", "tau", "rho",
          # parameters for transect random effects:
          "sigma_tr", "sigma_tr_sp", "delta_tr", "delta_tr_sp")

Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
rstan_options(auto_write = TRUE)

source_fn <- "SCAR_v3.stan"
output_fn <- "model_fit_v3.RData"

# Fit model --------------------------------------------------------------

smodel <- stan_model(source_fn, verbose = T)

fit <- sampling(
  smodel, data = stan_dat, pars = pars, refresh = 50,
  # iter = 50, warmup = 0, chains = 1, cores = 1, # to test before long run
  iter = 8000, warmup = 2000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.97, max_treedepth = 10)
)

save(fit, file = output_fn)
