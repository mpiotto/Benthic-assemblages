/*
  In this code we use the non-centered parameterization of random effects.
  Instead of defining
    betas ~ multi_normal(m, S),
  we define the parameter vector betas_raw:
    betas_raw ~ std_normal() // std_normal = normal(0, 1)
    betas = m + S_chol * betas_raw (S_chol es el factor de Cholesky)
    S_chol = cholesky_decompose(S)

  This is the multivariate analog of sampling from a normal distribution
  defined as
    x ~ normal(mu, sigma)
  using
    x_raw ~ normal(0, 1) // same as ~ std_normal()
    x = mu + sigma * x_raw.

  In addition, we never compute the S_chol matrix, which is expensive, but
  insted compute the (S_chol * betas_raw) product directly using kron_mult().
*/

functions {
  // Efficient Kronecker product matrix-vector multiplication, for square
  // matrices A and B, with the same dimension as x.
  vector kron_mult(matrix A, matrix B, vector x) {
    int rA = rows(A);
    int cA = cols(A);
    int rB = rows(B);
    int cB = cols(B);

    int n = rA * rB; // Total output size
    vector[n] y;  // Resulting vector

    matrix[rB, cA] X;  // Temporary matrix

    // Reshape x into an rB x cA matrix
    for (i in 1:cA) {
      X[, i] = segment(x, (i - 1) * rB + 1, rB);
    }

    // Compute the matrix product B * X
    X = B * X;

    // Flatten the (B * X) * A^T into a vector
    for (i in 1:rA) {
      y[((i - 1) * rB + 1):(i * rB)] = X * A[i]';
    }

    return y;
  }

  // Kronecker product (not used)
  matrix kronecker(matrix A, matrix B) {
    matrix[rows(A)*rows(B), cols(A)*cols(B)] kron;
    for (i in 1:cols(A)) {
      for (j in 1:rows(A)) {
        kron[((j-1)*rows(B)+1):(j*rows(B)), ((i-1)*cols(B)+1):(i*cols(B))] = A[j,i] * B;
      }
    }
    return kron;
  }
}

data {
  int<lower=1> N;            // number of photos
  int<lower=1> Nt;           // number of transects (random intercept)
  int<lower=1> Kpres;        // number of predictors for presence (int, y98, y09, y21, depth, s98, s09, s22)
  int<lower=1> Kcov;         // number of predictors for cover (int, y98, y09, y21, depth, s98, s09, s22)
  matrix[N, Kpres] Xpres;    // obs-level design matrix for environmental covariates
  matrix[N, Kcov] Xcov;      // obs-level design matrix for environmental covariates
  matrix[N, Nt] X_transect;  // obs-level design matrix for transect

  /*
    With this model structure 2 X matrices are not needed
    (same environmental covariates), but they are split in case that we want
    to include different covariates for presence and cover in future versions.
  */
  int<lower=1> N_J;     // num of genus
  int<lower=1> L_J;     // num group traits (habitat, rstrategy, fstrategy)
  matrix[N_J, L_J] TT;  // genus traits
  matrix[N_J, N_J] C;   // phylogenetic correlation matrix of genus

  // presence and cover
  int<lower=0, upper=1> Y_pres[N, N_J];  // binary presence matrix
  matrix[N, N_J] Y_cov;

  // data to evaluate cover likelihood efficiently
  int<lower=1> Npres; // total number of observations with presence
  int<lower=1> rows_cover[Npres]; // rows with cover (matrix indexation)
  int<lower=1> pres_length[N_J];  // number of rows with cover by genus
  int<lower=1> begin[N_J];

  // prior scale for z parameters
  real prior_sd_z;
}

transformed data {
  // Total amount of coefficients, to estimate the full vcov
  int<lower=1> K = Kpres + Kcov;

  // vector of 1s to compute phylogenetic correlation matrix
  vector[N_J] ones = rep_vector(1.0, N_J);

  // to compute the cover likelihood efficiently:
  int max_pres_length = max(pres_length);
}

/*
  The _raw parameters are those defined with normal(0, 1) prior and then scaled 
  (multip sigma) and translated (sum mu) as explained above. The transformation 
  is computed in the transformed parameters block.
*/

parameters {
  // Obtain genus-specific environmental responses

  // Correlation and scale to make variance-covariance matrices across parameter
  // types.
  vector<lower=0>[K] tau;
  cholesky_factor_corr[K] Omega_L;
  /*
     Correlation matrix for var-covar of betas (construct Sigma).
     We define the cholesky factor instead of Omega directly because it is more
     efficient and stable.
  */

  // coeffs for traits effects on species environmental responses (gamma)
  vector[L_J * K] z_raw;

  // genus environmental responses (for presence and abundance)
  vector[N_J * K] betas_raw;

  real<lower=0,upper=1> rho;  // phylogenetic signal of responses

  // Random effects of transect and transect * genus
  vector[Nt] delta_tr_raw;
  matrix[Nt, N_J] delta_tr_sp_raw;

  real<lower=0> sigma_tr;    // transects sd
  real<lower=0> sigma_tr_sp; // transects * genus sd

  // Cover precision parameter
  real<lower=0> phi;
}

