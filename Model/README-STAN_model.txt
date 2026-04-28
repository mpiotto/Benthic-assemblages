================================================================================
JOINT SPECIES DISTRIBUTION MODEL (JSDM) - Antarctic benthic assemblages
================================================================================

A hierarchical Joint Species Distribution Model (JSDM) implemented in Stan,
following the HMSC framework (Hierarchical Modelling of Species Communities).
It jointly models PRESENCE and COVER|PRESENCE of a set of Antarctic benthic species as
a function of environmental covariates, species functional traits and
phylogeny.


--------------------------------------------------------------------------------
1. OVERVIEW
--------------------------------------------------------------------------------

The model is a hurdle model with two components:

  1. PRESENCE - Bernoulli likelihood with logit link on the binary
     observations Y_pres.

  2. COVER - Beta likelihood (mean-precision parameterization) on the
     continuous observations Y_cov, evaluated ONLY in rows where presence
     is observed (Y_cov > 0).

Each species has its own vector of environmental responses, but these
responses are hierarchically structured by:

  - The species' functional traits (matrix TT), which determine the
    expected response.
  - The phylogeny among species (matrix C), which determines how similar
    the deviations from the expected value are.
  - A variance-covariance matrix among environmental coefficients, shared
    across species.


--------------------------------------------------------------------------------
2. MATHEMATICAL STRUCTURE
--------------------------------------------------------------------------------

2.1. Species-specific responses
-------------------------------

For each species j, the vector of environmental coefficients beta_j (of
length K = Kpres + Kcov) is modeled as:

    beta ~ MVN( m, S )
    m    = vec(TT * Z)
    S    = Sigma (kron) R

where:

  - Z       (L_J x K)     : effects of traits on environmental responses.
  - Sigma   (K x K)       : var-cov matrix among coefficients (shared
                            across species).
  - R       (N_J x N_J)   : phylogenetic correlation matrix:
                            R = rho * C + (1 - rho) * I.
  - rho     in [0, 1]     : phylogenetic signal of responses.


2.2. Linear predictor
---------------------

For each observation i and species j:

    logit(p_ij) = X_pres[i,] * b_pres[,j] + delta_tr[t(i)]
                                          + delta_tr_sp[t(i), j]
    logit(mu_ij) = X_cov[i,]  * b_cov[,j]

  - delta_tr     : random transect effect (shared across species).
  - delta_tr_sp  : random transect-by-species effect.
  - Random effects are applied ONLY in the presence component.


2.3. Likelihoods
----------------

    Y_pres[i,j]              ~ Bernoulli(p_ij)
    Y_cov[i,j] | Y_pres = 1  ~ Beta_proportion(mu_ij, phi)


--------------------------------------------------------------------------------
3. INPUT DATA (`data` block)
--------------------------------------------------------------------------------

  N             int                 Number of photos (observations).
  Nt            int                 Number of transects (random intercept).
  Kpres         int                 Number of predictors for presence (8).
  Kcov          int                 Number of predictors for cover (8).
  Xpres         matrix[N, Kpres]    Environmental design matrix (presence).
  Xcov          matrix[N, Kcov]     Environmental design matrix (cover).
  X_transect    matrix[N, Nt]       Design matrix for the transect effect.
  N_J           int                 Number of species.
  L_J           int                 Number of traits per species (habitat,
                                    r-strategy, f-strategy).
  TT            matrix[N_J, L_J]    Trait matrix.
  C             matrix[N_J, N_J]    Phylogenetic correlation matrix.
  Y_pres        int[N, N_J]         Binary presence matrix.
  Y_cov         matrix[N, N_J]      Cover matrix (proportion).
  Npres         int                 Total number of observations with
                                    presence.
  rows_cover    int[Npres]          Row indices with cover > 0 (concatenated
                                    by species).
  pres_length   int[N_J]            Number of rows with cover > 0 per
                                    species.
  begin         int[N_J]            Start index in rows_cover for each
                                    species.
  prior_sd_z    real                Prior scale for z.


Predictor structure (column order in `b`)
-----------------------------------------

Columns 1-8 correspond to presence and 9-16 to cover. For each block:

  Col   Meaning
  ----  ---------------------------------------------------------
  1     Intercept (reference: ~15 m depth, year 1994)
  2     Intercept differential, year 1998
  3     Intercept differential, year 2010
  4     Intercept differential, year 2022
  5     Slope with bathymetry
  6     Slope differential, year 1998
  7     Slope differential, year 2010
  8     Slope differential, year 2022


--------------------------------------------------------------------------------
4. ESTIMATED PARAMETERS
--------------------------------------------------------------------------------

  tau              vector[K]                       Scales of the var-cov
                                                   among coefficients.
  Omega_L          cholesky_factor_corr[K]         Cholesky factor of the
                                                   correlation among
                                                   coefficients.
  z_raw            vector[L_J * K]                 Trait effects on
                                                   responses (non-centered).
  betas_raw        vector[N_J * K]                 Species-specific
                                                   responses (non-centered).
  rho              real in [0, 1]                  Phylogenetic signal.
  delta_tr_raw     vector[Nt]                      Transect random effect
                                                   (non-centered).
  delta_tr_sp_raw  matrix[Nt, N_J]                 Transect-by-species
                                                   random effect
                                                   (non-centered).
  sigma_tr         real > 0                        SD of the transect
                                                   effect.
  sigma_tr_sp      real > 0                        SD of the
                                                   transect-by-species
                                                   effect.
  phi              real > 0                        Beta precision
                                                   parameter.


--------------------------------------------------------------------------------
5. PRIORS
--------------------------------------------------------------------------------

  Omega_L                          ~ lkj_corr_cholesky(2)
  tau                              ~ student_t(3, 0, 2.5)
  rho                              ~ beta(2, 2)
  z_raw, betas_raw                 ~ std_normal()
  delta_tr_raw, delta_tr_sp_raw    ~ std_normal()
  sigma_tr, sigma_tr_sp            ~ student_t(3, 0, 2)
  phi                              ~ inv_gamma(0.001, 0.001)


--------------------------------------------------------------------------------
6. COMPUTATIONAL TRICKS
--------------------------------------------------------------------------------

6.1. Non-centered parameterization
----------------------------------

All hierarchical effects are sampled in their non-centered form to improve
the geometry of the posterior and reduce HMC divergences. Instead of:

    beta ~ MVN(m, S)

the model uses:

    beta_raw ~ N(0, 1)
    beta     = m + chol(S) * beta_raw


6.2. Efficient Kronecker product (`kron_mult`)
----------------------------------------------

The matrix S = Sigma (kron) R would have dimension (N_J * K) x (N_J * K),
which is very expensive to build and decompose. The model leverages the
identity:

    (A (kron) B) * vec(X) = vec(B * X * A')

The function kron_mult(A, B, x) directly computes (A (kron) B) * x without
ever materializing the explicit Kronecker product. This drastically reduces
memory use and computation time.


6.3. Cover likelihood evaluated only on rows with presence
----------------------------------------------------------

Cover follows a Beta distribution and is only defined where presence is
observed. The model precomputes the indices rows_cover, pres_length and
begin in R/Python so that Stan, inside the species loop, evaluates the
likelihood ONLY on the relevant rows for each species, avoiding empty
calls and expensive masking.


--------------------------------------------------------------------------------
7. GENERATED QUANTITIES
--------------------------------------------------------------------------------

    Omega = Omega_L * Omega_L'

is the reconstructed correlation matrix among environmental coefficients
(useful for posterior interpretation).

--------------------------------------------------------------------------------
8. HOW TO RUN IT
--------------------------------------------------------------------------------

Minimal data structure (in R with `cmdstanr` or `rstan`):

    data_list <- list(
      N = nrow(Y_pres),
      Nt = nlevels(transect_id),
      Kpres = 8, Kcov = 8,
      Xpres = Xpres, Xcov = Xcov,
      X_transect = model.matrix(~ 0 + transect_id),
      N_J = ncol(Y_pres), L_J = ncol(TT),
      TT = TT, C = C,
      Y_pres = Y_pres, Y_cov = Y_cov,
      Npres = sum(Y_pres),
      rows_cover  = rows_cover,
      pres_length = pres_length,
      begin       = begin,
      prior_sd_z  = 1
    )

    fit <- model$sample(
      data = data_list,
      chains = 4, parallel_chains = 4,
      iter_warmup = 1000, iter_sampling = 1000,
      adapt_delta = 0.95
    )


Building rows_cover, pres_length and begin
------------------------------------------

    rows_cover  <- integer(0)
    pres_length <- integer(N_J)
    begin       <- integer(N_J)
    offset <- 1
    for (j in 1:N_J) {
      rj <- which(Y_pres[, j] == 1)
      rows_cover     <- c(rows_cover, rj)
      pres_length[j] <- length(rj)
      begin[j]       <- offset
      offset         <- offset + length(rj)
    }

