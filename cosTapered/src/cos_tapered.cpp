#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>
#include <BayesLogit.h>

#ifdef _OPENMP
# include <omp.h>
#endif

#include <climits>
#include <cmath>
#include <vector>

/*==========================================================================
  cos_tapered.cpp

  Functions:

    make_C_B_tapered()
      Direct/reference C_B^tap construction. This recomputes distances and
      taper weights every time. Keep this for testing and validation.

    prep_C_B_tapered_pairs()
      One-time precomputation for fixed H_BA, coordinates, gamma, and taper.
      This is analogous to building a sparse matrix pattern. It returns the
      contributing plot-cell pairs with d < gamma.

    make_C_B_tapered_from_pairs()
      Fast C_B^tap construction from precomputed pairs. Use this inside MCMC.

  Taper codes:
    1 = Wendland:  (1-r)^4 (1+4r) I(r < 1)
    2 = spherical: 1 - 1.5r + 0.5r^3 I(r < 1)

  Compile:
    R CMD SHLIB cos_tapered.cpp
==========================================================================*/

static double as_real_scalar(SEXP x, const char *where)
{
  if(x == R_NilValue || LENGTH(x) != 1)
    Rf_error("%s must be a scalar", where);

  if(TYPEOF(x) == REALSXP)
    return REAL(x)[0];

  if(TYPEOF(x) == INTSXP)
    return (double)INTEGER(x)[0];

  Rf_error("%s must be numeric or integer", where);
  return 0.0;
}

static int as_int_scalar(SEXP x, const char *where)
{
  if(x == R_NilValue || LENGTH(x) != 1)
    Rf_error("%s must be a scalar", where);

  if(TYPEOF(x) == INTSXP)
    return INTEGER(x)[0];

  if(TYPEOF(x) == REALSXP)
    return (int)REAL(x)[0];

  Rf_error("%s must be integer or numeric", where);
  return 0;
}

static void require_real_vector(SEXP x, const char *where)
{
  if(TYPEOF(x) != REALSXP)
    Rf_error("%s must be numeric", where);
}

static void require_integer_vector(SEXP x, const char *where)
{
  if(TYPEOF(x) != INTSXP)
    Rf_error("%s must be integer", where);
}

static void check_common_inputs(SEXP plot_start_r,
                                SEXP h_r,
                                SEXP x_r,
                                SEXP y_r,
                                double phi,
                                double gamma,
                                int taper_code,
                                int *n_b_out,
                                int *n_nz_out)
{
  require_integer_vector(plot_start_r, "plot_start");
  require_real_vector(h_r, "h");
  require_real_vector(x_r, "x");
  require_real_vector(y_r, "y");

  if(!R_finite(phi) || phi <= 0.0)
    Rf_error("phi must be finite and positive");
  if(!R_finite(gamma) || gamma <= 0.0)
    Rf_error("gamma must be finite and positive");
  if(taper_code != 1 && taper_code != 2)
    Rf_error("taper_code must be 1 (Wendland) or 2 (spherical)");

  int n_b = LENGTH(plot_start_r) - 1;
  int n_nz = LENGTH(h_r);

  if(n_b <= 0)
    Rf_error("plot_start must have length at least 2");
  if(LENGTH(x_r) != n_nz || LENGTH(y_r) != n_nz)
    Rf_error("h, x, and y must have the same length");

  int *plot_start = INTEGER(plot_start_r);

  if(plot_start[0] != 1)
    Rf_error("plot_start must use R-style indexing and start at 1");
  if(plot_start[n_b] != n_nz + 1)
    Rf_error("last plot_start must equal length(h) + 1");

  for(int l = 0; l < n_b; l++){
    if(plot_start[l] < 1 || plot_start[l + 1] < plot_start[l])
      Rf_error("plot_start must be nondecreasing and positive");
  }

  *n_b_out = n_b;
  *n_nz_out = n_nz;
}

static double wendland_taper(double d, double gamma)
{
  if(d >= gamma)
    return 0.0;

  double r = d / gamma;
  double omr = 1.0 - r;

  return std::pow(omr, 4.0) * (1.0 + 4.0 * r);
}

static double spherical_taper(double d, double gamma)
{
  if(d >= gamma)
    return 0.0;

  double r = d / gamma;

  return 1.0 - 1.5 * r + 0.5 * r * r * r;
}

static double taper_value(double d, double gamma, int taper_code)
{
  if(taper_code == 1)
    return wendland_taper(d, gamma);

  return spherical_taper(d, gamma);
}

static void require_real_matrix(SEXP x, const char *where)
{
  require_real_vector(x, where);

  SEXP dim_r = Rf_getAttrib(x, R_DimSymbol);
  if(dim_r == R_NilValue || LENGTH(dim_r) != 2)
    Rf_error("%s must be a numeric matrix", where);
}

static double phi_from_z(double z, double phi_lower, double phi_upper)
{
  return phi_upper - (phi_upper - phi_lower) / (1.0 + std::exp(z));
}

static void fill_C_B_from_pairs_cpp(double *C_B,
                                    int n_b,
                                    const int *pair_l,
                                    const int *pair_k,
                                    const double *pair_d,
                                    const double *pair_wtap,
                                    R_xlen_t n_pair,
                                    double phi,
                                    int n_threads)
{
  R_xlen_t n_elem = static_cast<R_xlen_t>(n_b) * static_cast<R_xlen_t>(n_b);

  for(R_xlen_t i = 0; i < n_elem; i++)
    C_B[i] = 0.0;

#ifdef _OPENMP
  if(n_threads > 1)
    omp_set_num_threads(n_threads);

#pragma omp parallel if(n_threads > 1)
  {
    std::vector<double> C_local(n_elem, 0.0);

#pragma omp for schedule(static)
    for(R_xlen_t r = 0; r < n_pair; r++){
      int l = pair_l[r] - 1;
      int k = pair_k[r] - 1;

      if(l < 0 || l >= n_b || k < 0 || k >= n_b)
        continue;

      double val = pair_wtap[r] * std::exp(-phi * pair_d[r]);

      C_local[l + n_b * k] += val;

      if(l != k)
        C_local[k + n_b * l] += val;
    }

#pragma omp critical
    {
      for(R_xlen_t i = 0; i < n_elem; i++)
        C_B[i] += C_local[i];
    }
  }
#else
  for(R_xlen_t r = 0; r < n_pair; r++){
    int l = pair_l[r] - 1;
    int k = pair_k[r] - 1;

    if(l < 0 || l >= n_b || k < 0 || k >= n_b)
      Rf_error("pair_l or pair_k contains an out-of-range plot index");

    double val = pair_wtap[r] * std::exp(-phi * pair_d[r]);

    C_B[l + n_b * k] += val;

    if(l != k)
      C_B[k + n_b * l] += val;
  }
#endif
}

static bool chol_upper_in_place(std::vector<double> &A, int n)
{
  char uplo = 'U';
  int info = 0;

  F77_CALL(dpotrf)(&uplo, &n, A.data(), &n, &info FCONE);

  if(info != 0)
    return false;

  for(int j = 0; j < n; j++){
    for(int i = j + 1; i < n; i++)
      A[i + n * j] = 0.0;
  }

  return true;
}

static bool chol_inverse_upper_in_place(std::vector<double> &A, int n)
{
  if(!chol_upper_in_place(A, n))
    return false;

  char uplo = 'U';
  int info = 0;
  F77_CALL(dpotri)(&uplo, &n, A.data(), &n, &info FCONE);

  if(info != 0)
    return false;

  for(int j = 0; j < n; j++){
    for(int i = j + 1; i < n; i++)
      A[i + n * j] = A[j + n * i];
  }

  return true;
}

static void solve_chol_upper_in_place(const std::vector<double> &R,
                                      int n,
                                      std::vector<double> &b)
{
  char uplo = 'U';
  int nrhs = 1;
  int info = 0;

  F77_CALL(dpotrs)(&uplo, &n, &nrhs, R.data(), &n, b.data(), &n, &info FCONE);

  if(info != 0)
    Rf_error("Cholesky solve failed");
}

static void rmvn_from_precision(const std::vector<double> &Q,
                                const std::vector<double> &b,
                                int d,
                                std::vector<double> &draw,
                                std::vector<double> &mean)
{
  std::vector<double> R = Q;
  if(!chol_upper_in_place(R, d))
    Rf_error("Conditional precision is not positive definite");

  mean = b;
  solve_chol_upper_in_place(R, d, mean);

  std::vector<double> z(d, 0.0);
  for(int i = 0; i < d; i++)
    z[i] = norm_rand();

  draw = mean;
  char trans = 'N';
  char diag = 'N';
  char uplo = 'U';
  int inc = 1;
  F77_CALL(dtrsv)(&uplo, &trans, &diag, &d, R.data(), &d, z.data(), &inc FCONE FCONE FCONE);

  for(int i = 0; i < d; i++)
    draw[i] += z[i];
}

static double log_target_gaussian_cpp(const double *theta,
                                      int d,
                                      bool spatial,
                                      const double *y_B,
                                      const double *mu_y,
                                      const double *XVX,
                                      const double *D_h,
                                      int n_b,
                                      const int *pair_l,
                                      const int *pair_k,
                                      const double *pair_d,
                                      const double *pair_wtap,
                                      R_xlen_t n_pair,
                                      double d_h_bar,
                                      double tau_B_shape,
                                      double tau_B_scale,
                                      double sigma_shape,
                                      double sigma_scale,
                                      double phi_lower,
                                      double phi_upper,
                                      int n_threads,
                                      std::vector<double> &V,
                                      std::vector<double> &C_B,
                                      std::vector<double> &r)
{
  double log_tau_B_sq = theta[0];
  double tau_B_sq = std::exp(log_tau_B_sq);

  if(!R_finite(tau_B_sq) || tau_B_sq <= 0.0)
    return R_NegInf;

  double tau_sq = tau_B_sq / d_h_bar;
  double sigma_sq = R_NaReal;
  double phi = R_NaReal;

  if(spatial){
    if(d < 3)
      return R_NegInf;

    sigma_sq = std::exp(theta[1]);
    phi = phi_from_z(theta[2], phi_lower, phi_upper);

    if(!R_finite(sigma_sq) || sigma_sq <= 0.0 ||
       !R_finite(phi) || phi <= phi_lower || phi >= phi_upper)
      return R_NegInf;
  }

  R_xlen_t n_elem = static_cast<R_xlen_t>(n_b) * static_cast<R_xlen_t>(n_b);

  for(R_xlen_t i = 0; i < n_elem; i++)
    V[i] = XVX[i] + tau_sq * D_h[i];

  if(spatial){
    fill_C_B_from_pairs_cpp(C_B.data(), n_b, pair_l, pair_k, pair_d,
                            pair_wtap, n_pair, phi, n_threads);

    for(R_xlen_t i = 0; i < n_elem; i++)
      V[i] += sigma_sq * C_B[i];
  }

  if(!chol_upper_in_place(V, n_b))
    return R_NegInf;

  for(int i = 0; i < n_b; i++)
    r[i] = y_B[i] - mu_y[i];

  char uplo = 'U';
  int nrhs = 1;
  int info = 0;
  F77_CALL(dpotrs)(&uplo, &n_b, &nrhs, V.data(), &n_b, r.data(), &n_b, &info FCONE);

  if(info != 0)
    return R_NegInf;

  double quad = 0.0;
  for(int i = 0; i < n_b; i++)
    quad += (y_B[i] - mu_y[i]) * r[i];

  double log_det = 0.0;
  for(int i = 0; i < n_b; i++){
    double rii = V[i + n_b * i];
    if(!R_finite(rii) || rii <= 0.0)
      return R_NegInf;
    log_det += 2.0 * std::log(rii);
  }

  double log_lik = -0.5 * (
    static_cast<double>(n_b) * std::log(2.0 * M_PI) + log_det + quad
  );

  double log_prior_tau_B =
    tau_B_shape * std::log(tau_B_scale) -
    lgammafn(tau_B_shape) -
    (tau_B_shape + 1.0) * std::log(tau_B_sq) -
    tau_B_scale / tau_B_sq;

  double out = log_lik + log_prior_tau_B + std::log(tau_B_sq);

  if(spatial){
    double log_prior_sigma =
      sigma_shape * std::log(sigma_scale) -
      lgammafn(sigma_shape) -
      (sigma_shape + 1.0) * std::log(sigma_sq) -
      sigma_scale / sigma_sq;

    double log_prior_phi = -std::log(phi_upper - phi_lower);
    double log_jac_extra =
      std::log(sigma_sq) +
      std::log(phi - phi_lower) +
      std::log(phi_upper - phi) -
      std::log(phi_upper - phi_lower);

    out += log_prior_sigma + log_prior_phi + log_jac_extra;
  }

  if(!R_finite(out))
    return R_NegInf;

  return out;
}

static void covariance_from_batch(const std::vector<double> &batch_samples,
                                  int batch_length,
                                  int d,
                                  std::vector<double> &Sigma_hat)
{
  std::vector<double> means(d, 0.0);

  for(int j = 0; j < batch_length; j++){
    for(int k = 0; k < d; k++)
      means[k] += batch_samples[j + batch_length * k];
  }

  for(int k = 0; k < d; k++)
    means[k] /= static_cast<double>(batch_length);

  for(int i = 0; i < d * d; i++)
    Sigma_hat[i] = 0.0;

  double denom = static_cast<double>(batch_length - 1);
  for(int j = 0; j < batch_length; j++){
    for(int a = 0; a < d; a++){
      double da = batch_samples[j + batch_length * a] - means[a];
      for(int b = 0; b <= a; b++){
        double db = batch_samples[j + batch_length * b] - means[b];
        Sigma_hat[a + d * b] += da * db / denom;
      }
    }
  }

  for(int a = 0; a < d; a++){
    for(int b = 0; b < a; b++)
      Sigma_hat[b + d * a] = Sigma_hat[a + d * b];
  }
}

/* Adaptive random-walk Metropolis for the current Gaussian COS target. */

extern "C" SEXP gaussian_metrop_sampler(SEXP y_B_r,
                                         SEXP mu_y_r,
                                         SEXP XVX_r,
                                         SEXP D_h_r,
                                         SEXP pair_l_r,
                                         SEXP pair_k_r,
                                         SEXP pair_d_r,
                                         SEXP pair_wtap_r,
                                         SEXP n_b_pairs_r,
                                         SEXP spatial_r,
                                         SEXP starting_r,
                                         SEXP tuning_r,
                                         SEXP n_batch_r,
                                         SEXP batch_length_r,
                                         SEXP accept_rate_r,
                                         SEXP c0_r,
                                         SEXP c1_r,
                                         SEXP report_r,
                                         SEXP verbose_r,
                                         SEXP d_h_bar_r,
                                         SEXP tau_B_shape_r,
                                         SEXP tau_B_scale_r,
                                         SEXP sigma_shape_r,
                                         SEXP sigma_scale_r,
                                         SEXP phi_lower_r,
                                         SEXP phi_upper_r,
                                         SEXP n_threads_r)
{
  require_real_vector(y_B_r, "y_B");
  require_real_vector(mu_y_r, "mu_y");
  require_real_matrix(XVX_r, "X_B_V_beta_X_B");
  require_real_matrix(D_h_r, "D_h");
  require_real_vector(starting_r, "starting");
  require_real_vector(tuning_r, "tuning");

  int n_b = LENGTH(y_B_r);
  int d = LENGTH(starting_r);

  if(LENGTH(mu_y_r) != n_b)
    Rf_error("mu_y must have the same length as y_B");
  if(LENGTH(tuning_r) != d)
    Rf_error("tuning must have the same length as starting");

  SEXP dim_xvx = Rf_getAttrib(XVX_r, R_DimSymbol);
  SEXP dim_dh = Rf_getAttrib(D_h_r, R_DimSymbol);

  if(INTEGER(dim_xvx)[0] != n_b || INTEGER(dim_xvx)[1] != n_b)
    Rf_error("X_B_V_beta_X_B must be n_B by n_B");
  if(INTEGER(dim_dh)[0] != n_b || INTEGER(dim_dh)[1] != n_b)
    Rf_error("D_h must be n_B by n_B");

  bool spatial = Rf_asLogical(spatial_r) == TRUE;

  const int *pair_l = NULL;
  const int *pair_k = NULL;
  const double *pair_d = NULL;
  const double *pair_wtap = NULL;
  R_xlen_t n_pair = 0;

  if(spatial){
    require_integer_vector(pair_l_r, "pair_l");
    require_integer_vector(pair_k_r, "pair_k");
    require_real_vector(pair_d_r, "pair_d");
    require_real_vector(pair_wtap_r, "pair_wtap");

    n_pair = XLENGTH(pair_l_r);
    if(XLENGTH(pair_k_r) != n_pair ||
       XLENGTH(pair_d_r) != n_pair ||
       XLENGTH(pair_wtap_r) != n_pair)
      Rf_error("pair_l, pair_k, pair_d, and pair_wtap must have the same length");

    int n_b_pairs = as_int_scalar(n_b_pairs_r, "n_b_pairs");
    if(n_b_pairs != n_b)
      Rf_error("C_B_pairs$n_b must match length(y_B)");

    pair_l = INTEGER(pair_l_r);
    pair_k = INTEGER(pair_k_r);
    pair_d = REAL(pair_d_r);
    pair_wtap = REAL(pair_wtap_r);

    for(R_xlen_t r = 0; r < n_pair; r++){
      if(pair_l[r] < 1 || pair_l[r] > n_b ||
         pair_k[r] < 1 || pair_k[r] > n_b)
        Rf_error("pair_l or pair_k contains an out-of-range plot index");
      if(!R_finite(pair_d[r]) || pair_d[r] < 0.0)
        Rf_error("pair_d must contain finite nonnegative distances");
      if(!R_finite(pair_wtap[r]))
        Rf_error("pair_wtap must contain finite values");
    }
  }

  int n_batch = as_int_scalar(n_batch_r, "n_batch");
  int batch_length = as_int_scalar(batch_length_r, "batch_length");
  int report = as_int_scalar(report_r, "report");
  int n_threads = as_int_scalar(n_threads_r, "n_threads");
  bool verbose = Rf_asLogical(verbose_r) == TRUE;

  double accept_rate = as_real_scalar(accept_rate_r, "accept_rate");
  double c0 = as_real_scalar(c0_r, "c0");
  double c1 = as_real_scalar(c1_r, "c1");
  double d_h_bar = as_real_scalar(d_h_bar_r, "d_h_bar");
  double tau_B_shape = as_real_scalar(tau_B_shape_r, "tau_B_shape");
  double tau_B_scale = as_real_scalar(tau_B_scale_r, "tau_B_scale");
  double sigma_shape = as_real_scalar(sigma_shape_r, "sigma_shape");
  double sigma_scale = as_real_scalar(sigma_scale_r, "sigma_scale");
  double phi_lower = as_real_scalar(phi_lower_r, "phi_lower");
  double phi_upper = as_real_scalar(phi_upper_r, "phi_upper");

  if(n_batch < 1 || batch_length < 1)
    Rf_error("n_batch and batch_length must be positive");
  if(report < 1)
    report = n_batch + 1;
  if(d < 1)
    Rf_error("starting must contain at least one value");
  if(spatial && d != 3)
    Rf_error("spatial Gaussian sampler expects three parameters");
  if(!spatial && d != 1)
    Rf_error("nonspatial Gaussian sampler expects one parameter");
  if(n_threads < 1)
    n_threads = 1;

  int n_save = n_batch * batch_length;

  SEXP samples_r, lp_samples_r, accept_r, batch_accept_rate_r;
  SEXP sigma_sq_m_r, Sigma0_r, proposal_cov_r, out_r, names_r;

  PROTECT(samples_r = Rf_allocMatrix(REALSXP, n_save, d));
  PROTECT(lp_samples_r = Rf_allocVector(REALSXP, n_save));
  PROTECT(accept_r = Rf_allocMatrix(INTSXP, n_batch, batch_length));
  PROTECT(batch_accept_rate_r = Rf_allocVector(REALSXP, n_batch));
  PROTECT(sigma_sq_m_r = Rf_allocVector(REALSXP, 1));
  PROTECT(Sigma0_r = Rf_allocMatrix(REALSXP, d, d));
  PROTECT(proposal_cov_r = Rf_allocMatrix(REALSXP, d, d));

  double *samples = REAL(samples_r);
  double *lp_samples = REAL(lp_samples_r);
  int *accept = INTEGER(accept_r);
  double *batch_accept_rate = REAL(batch_accept_rate_r);
  double *Sigma0_out = REAL(Sigma0_r);
  double *proposal_cov = REAL(proposal_cov_r);

  for(int i = 0; i < n_save * d; i++)
    samples[i] = NA_REAL;
  for(int i = 0; i < n_save; i++)
    lp_samples[i] = NA_REAL;
  for(int i = 0; i < n_batch * batch_length; i++)
    accept[i] = 0;

  std::vector<double> theta(d);
  std::vector<double> theta_prop(d);
  std::vector<double> tuning(d);
  std::vector<double> Sigma0(d * d, 0.0);
  std::vector<double> Sigma_hat(d * d, 0.0);
  std::vector<double> R_prop(d * d, 0.0);
  std::vector<double> z(d, 0.0);
  std::vector<double> V(static_cast<R_xlen_t>(n_b) * n_b, 0.0);
  std::vector<double> C_B(static_cast<R_xlen_t>(n_b) * n_b, 0.0);
  std::vector<double> r(n_b, 0.0);
  std::vector<double> batch_samples(batch_length * d, 0.0);

  for(int i = 0; i < d; i++){
    theta[i] = REAL(starting_r)[i];
    tuning[i] = REAL(tuning_r)[i];
    Sigma0[i + d * i] = tuning[i] * tuning[i];
  }

  double sigma_sq_m = 2.4 * 2.4 / static_cast<double>(d);
  int save_i = 0;

  double lp = log_target_gaussian_cpp(
    theta.data(), d, spatial, REAL(y_B_r), REAL(mu_y_r), REAL(XVX_r),
    REAL(D_h_r), n_b, pair_l, pair_k, pair_d, pair_wtap, n_pair,
    d_h_bar, tau_B_shape, tau_B_scale, sigma_shape, sigma_scale,
    phi_lower, phi_upper, n_threads, V, C_B, r
  );

  if(!R_finite(lp))
    Rf_error("Starting value has non-finite log target density");

  GetRNGstate();

  for(int b = 0; b < n_batch; b++){
    for(int i = 0; i < d * d; i++)
      R_prop[i] = sigma_sq_m * Sigma0[i];

    if(!chol_upper_in_place(R_prop, d)){
      PutRNGstate();
      Rf_error("Proposal covariance is not positive definite");
    }

    int batch_accept = 0;

    for(int j = 0; j < batch_length; j++){
      for(int k = 0; k < d; k++)
        z[k] = norm_rand();

      for(int k = 0; k < d; k++){
        double step = 0.0;
        for(int a = 0; a <= k; a++)
          step += R_prop[a + d * k] * z[a];
        theta_prop[k] = theta[k] + step;
      }

      double lp_prop = log_target_gaussian_cpp(
        theta_prop.data(), d, spatial, REAL(y_B_r), REAL(mu_y_r),
        REAL(XVX_r), REAL(D_h_r), n_b, pair_l, pair_k, pair_d,
        pair_wtap, n_pair, d_h_bar, tau_B_shape, tau_B_scale,
        sigma_shape, sigma_scale, phi_lower, phi_upper, n_threads,
        V, C_B, r
      );

      if(std::log(unif_rand()) < lp_prop - lp){
        for(int k = 0; k < d; k++)
          theta[k] = theta_prop[k];
        lp = lp_prop;
        accept[b + n_batch * j] = 1;
        batch_accept++;
      }

      for(int k = 0; k < d; k++){
        samples[save_i + n_save * k] = theta[k];
        batch_samples[j + batch_length * k] = theta[k];
      }
      lp_samples[save_i] = lp;
      save_i++;
    }

    double rhat = static_cast<double>(batch_accept) /
      static_cast<double>(batch_length);
    batch_accept_rate[b] = rhat;

    double gamma1 = 1.0 / std::pow(static_cast<double>(b + 2), c1);
    double gamma2 = c0 * gamma1;

    sigma_sq_m = std::exp(std::log(sigma_sq_m) +
      gamma2 * (rhat - accept_rate));

    if(batch_length > 1){
      covariance_from_batch(batch_samples, batch_length, d, Sigma_hat);

      bool ok = true;
      for(int i = 0; i < d * d; i++){
        if(!R_finite(Sigma_hat[i])){
          ok = false;
          break;
        }
      }

      if(ok){
        for(int i = 0; i < d * d; i++)
          Sigma0[i] += gamma1 * (Sigma_hat[i] - Sigma0[i]);
      }
    }

    if(verbose && (b == 0 || ((b + 1) % report == 0))){
      Rprintf("batch %d accept = %.3f sigma_sq_m = %.4g lp = %.6g\n",
              b + 1, rhat, sigma_sq_m, lp);
    }

    if((b + 1) % 10 == 0)
      R_CheckUserInterrupt();
  }

  PutRNGstate();

  REAL(sigma_sq_m_r)[0] = sigma_sq_m;

  for(int i = 0; i < d * d; i++){
    Sigma0_out[i] = Sigma0[i];
    proposal_cov[i] = sigma_sq_m * Sigma0[i];
  }

  PROTECT(out_r = Rf_allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out_r, 0, samples_r);
  SET_VECTOR_ELT(out_r, 1, lp_samples_r);
  SET_VECTOR_ELT(out_r, 2, accept_r);
  SET_VECTOR_ELT(out_r, 3, batch_accept_rate_r);
  SET_VECTOR_ELT(out_r, 4, sigma_sq_m_r);
  SET_VECTOR_ELT(out_r, 5, Sigma0_r);
  SET_VECTOR_ELT(out_r, 6, proposal_cov_r);

  PROTECT(names_r = Rf_allocVector(STRSXP, 7));
  SET_STRING_ELT(names_r, 0, Rf_mkChar("p.theta.samples"));
  SET_STRING_ELT(names_r, 1, Rf_mkChar("p.lp.samples"));
  SET_STRING_ELT(names_r, 2, Rf_mkChar("accept"));
  SET_STRING_ELT(names_r, 3, Rf_mkChar("batch.accept.rate"));
  SET_STRING_ELT(names_r, 4, Rf_mkChar("sigma_sq_m"));
  SET_STRING_ELT(names_r, 5, Rf_mkChar("Sigma0"));
  SET_STRING_ELT(names_r, 6, Rf_mkChar("proposal.cov"));

  Rf_setAttrib(out_r, R_NamesSymbol, names_r);

  UNPROTECT(9);
  return out_r;
}

static double log_theta_binomial_spatial(const double *theta,
                                         const double *omega_B,
                                         int n_b,
                                         const int *pair_l,
                                         const int *pair_k,
                                         const double *pair_d,
                                         const double *pair_wtap,
                                         R_xlen_t n_pair,
                                         double sigma_shape,
                                         double sigma_scale,
                                         double phi_lower,
                                         double phi_upper,
                                         int n_threads,
                                         std::vector<double> &C_B,
                                         std::vector<double> &v)
{
  double sigma_sq = std::exp(theta[0]);
  double phi = phi_from_z(theta[1], phi_lower, phi_upper);

  if(!R_finite(sigma_sq) || sigma_sq <= 0.0 ||
     !R_finite(phi) || phi <= phi_lower || phi >= phi_upper)
    return R_NegInf;

  fill_C_B_from_pairs_cpp(C_B.data(), n_b, pair_l, pair_k, pair_d,
                          pair_wtap, n_pair, phi, n_threads);

  if(!chol_upper_in_place(C_B, n_b))
    return R_NegInf;

  for(int i = 0; i < n_b; i++)
    v[i] = omega_B[i];

  solve_chol_upper_in_place(C_B, n_b, v);

  double quad = 0.0;
  for(int i = 0; i < n_b; i++)
    quad += omega_B[i] * v[i];

  double log_det_C = 0.0;
  for(int i = 0; i < n_b; i++)
    log_det_C += 2.0 * std::log(C_B[i + n_b * i]);

  double log_prior_omega = -0.5 * (
    static_cast<double>(n_b) * std::log(2.0 * M_PI * sigma_sq) +
      log_det_C + quad / sigma_sq
  );

  double log_prior_sigma =
    sigma_shape * std::log(sigma_scale) -
    lgammafn(sigma_shape) -
    (sigma_shape + 1.0) * std::log(sigma_sq) -
    sigma_scale / sigma_sq;

  double log_prior_phi = -std::log(phi_upper - phi_lower);
  double log_jac =
    std::log(sigma_sq) +
    std::log(phi - phi_lower) +
    std::log(phi_upper - phi) -
    std::log(phi_upper - phi_lower);

  double out = log_prior_omega + log_prior_sigma + log_prior_phi + log_jac;
  if(!R_finite(out))
    return R_NegInf;

  return out;
}

static void update_pg_response(const double *y,
                               const int *trials,
                               const double *offset,
                               int family_code,
                               double nb_size,
                               const double *X,
                               int n_b,
                               int p,
                               bool spatial,
                               const std::vector<double> &beta,
                               const std::vector<double> &omega_B,
                               BayesLogit_rpg_hybrid_t pg,
                               std::vector<double> &pg_w,
                               std::vector<double> &z)
{
  for(int i = 0; i < n_b; i++){
    double eta = offset[i];
    for(int j = 0; j < p; j++)
      eta += X[i + n_b * j] * beta[j];
    if(spatial)
      eta += omega_B[i];

    double shape = 0.0;
    double kappa = 0.0;
    double shift = 0.0;
    double tilt = eta;

    if(family_code == 1){
      shape = static_cast<double>(trials[i]);
      kappa = y[i] - 0.5 * shape;
    } else if(family_code == 2){
      shape = y[i] + nb_size;
      kappa = 0.5 * (y[i] - nb_size);
      shift = std::log(nb_size);
      tilt = eta - shift;
    } else {
      Rf_error("Unsupported Polya-Gamma family code");
    }

    double omega = pg(shape, tilt);

    if(!R_finite(omega) || omega <= 0.0)
      Rf_error("Polya-Gamma draw produced a non-positive observation precision");

    pg_w[i] = omega;
    z[i] = kappa / omega + shift - offset[i];
  }
}

static void sample_pg_alpha(const double *X,
                            int n_b,
                            int p,
                            bool spatial,
                            const std::vector<double> &pg_w,
                            const std::vector<double> &z,
                            const double *mu_beta,
                            const double *V_beta_inv,
                            const std::vector<double> &C_B_inv,
                            double sigma_sq,
                            std::vector<double> &beta,
                            std::vector<double> &omega_B,
                            std::vector<double> &alpha_mean)
{
  int d = p + (spatial ? n_b : 0);
  std::vector<double> Q(static_cast<R_xlen_t>(d) * d, 0.0);
  std::vector<double> b(d, 0.0);
  std::vector<double> alpha_draw(d, 0.0);

  for(int i = 0; i < n_b; i++){
    double wi = pg_w[i];
    double wz = wi * z[i];

    for(int a = 0; a < p; a++){
      double xia = X[i + n_b * a];
      b[a] += xia * wz;

      for(int c = 0; c < p; c++)
        Q[a + d * c] += xia * wi * X[i + n_b * c];

      if(spatial){
        int oi = p + i;
        Q[a + d * oi] += xia * wi;
        Q[oi + d * a] += xia * wi;
      }
    }

    if(spatial){
      int oi = p + i;
      b[oi] += wz;
      Q[oi + d * oi] += wi;
    }
  }

  for(int a = 0; a < p; a++){
    for(int c = 0; c < p; c++)
      Q[a + d * c] += V_beta_inv[a + p * c];

    for(int c = 0; c < p; c++)
      b[a] += V_beta_inv[a + p * c] * mu_beta[c];
  }

  if(spatial){
    for(int a = 0; a < n_b; a++){
      for(int c = 0; c < n_b; c++)
        Q[(p + a) + d * (p + c)] += C_B_inv[a + n_b * c] / sigma_sq;
    }
  }

  rmvn_from_precision(Q, b, d, alpha_draw, alpha_mean);

  for(int j = 0; j < p; j++)
    beta[j] = alpha_draw[j];

  if(spatial){
    for(int i = 0; i < n_b; i++)
      omega_B[i] = alpha_draw[p + i];
  } else {
    for(int i = 0; i < n_b; i++)
      omega_B[i] = 0.0;
  }
}

extern "C" SEXP binomial_pg_sampler(SEXP y_B_r,
                                     SEXP trials_r,
                                     SEXP offset_r,
                                     SEXP family_code_r,
                                     SEXP nb_size_r,
                                     SEXP X_B_r,
                                     SEXP pair_l_r,
                                     SEXP pair_k_r,
                                     SEXP pair_d_r,
                                     SEXP pair_wtap_r,
                                     SEXP n_b_pairs_r,
                                     SEXP spatial_r,
                                     SEXP beta_mu_r,
                                     SEXP V_beta_r,
                                     SEXP starting_beta_r,
                                     SEXP starting_sigma_sq_r,
                                     SEXP starting_phi_r,
                                     SEXP tuning_r,
                                     SEXP n_batch_r,
                                     SEXP batch_length_r,
                                     SEXP accept_rate_r,
                                     SEXP c0_r,
                                     SEXP c1_r,
                                     SEXP report_r,
                                     SEXP verbose_r,
                                     SEXP sigma_shape_r,
                                     SEXP sigma_scale_r,
                                     SEXP phi_lower_r,
                                     SEXP phi_upper_r,
                                     SEXP n_threads_r)
{
  require_real_vector(y_B_r, "y_B");
  require_integer_vector(trials_r, "trials");
  require_real_vector(offset_r, "offset");
  require_real_matrix(X_B_r, "X_B");
  require_real_vector(beta_mu_r, "beta_mu");
  require_real_matrix(V_beta_r, "V_beta");
  require_real_vector(starting_beta_r, "starting_beta");
  require_real_vector(tuning_r, "tuning");

  int n_b = LENGTH(y_B_r);
  int p = LENGTH(beta_mu_r);
  bool spatial = Rf_asLogical(spatial_r) == TRUE;
  int family_code = as_int_scalar(family_code_r, "family_code");
  double nb_size = as_real_scalar(nb_size_r, "nb_size");

  SEXP dim_x = Rf_getAttrib(X_B_r, R_DimSymbol);
  SEXP dim_vb = Rf_getAttrib(V_beta_r, R_DimSymbol);
  if(INTEGER(dim_x)[0] != n_b || INTEGER(dim_x)[1] != p)
    Rf_error("X_B must be n_B by length(beta_mu)");
  if(INTEGER(dim_vb)[0] != p || INTEGER(dim_vb)[1] != p)
    Rf_error("V_beta must be p by p");
  if(LENGTH(trials_r) != n_b)
    Rf_error("trials must have length n_B");
  if(LENGTH(offset_r) != n_b)
    Rf_error("offset must have length n_B");
  if(LENGTH(starting_beta_r) != p)
    Rf_error("starting_beta must have length length(beta_mu)");
  if(family_code != 1 && family_code != 2)
    Rf_error("family_code must be 1 (binomial) or 2 (negative_binomial)");
  if(family_code == 2 && (!R_finite(nb_size) || nb_size <= 0.0))
    Rf_error("negative-binomial size must be positive and finite");

  const int *pair_l = NULL;
  const int *pair_k = NULL;
  const double *pair_d = NULL;
  const double *pair_wtap = NULL;
  R_xlen_t n_pair = 0;

  if(spatial){
    require_integer_vector(pair_l_r, "pair_l");
    require_integer_vector(pair_k_r, "pair_k");
    require_real_vector(pair_d_r, "pair_d");
    require_real_vector(pair_wtap_r, "pair_wtap");

    n_pair = XLENGTH(pair_l_r);
    if(XLENGTH(pair_k_r) != n_pair ||
       XLENGTH(pair_d_r) != n_pair ||
       XLENGTH(pair_wtap_r) != n_pair)
      Rf_error("pair_l, pair_k, pair_d, and pair_wtap must have the same length");

    int n_b_pairs = as_int_scalar(n_b_pairs_r, "n_b_pairs");
    if(n_b_pairs != n_b)
      Rf_error("C_B_pairs$n_b must match length(y_B)");

    pair_l = INTEGER(pair_l_r);
    pair_k = INTEGER(pair_k_r);
    pair_d = REAL(pair_d_r);
    pair_wtap = REAL(pair_wtap_r);

    for(R_xlen_t r = 0; r < n_pair; r++){
      if(pair_l[r] < 1 || pair_l[r] > n_b ||
         pair_k[r] < 1 || pair_k[r] > n_b)
        Rf_error("pair_l or pair_k contains an out-of-range plot index");
      if(!R_finite(pair_d[r]) || pair_d[r] < 0.0)
        Rf_error("pair_d must contain finite nonnegative distances");
      if(!R_finite(pair_wtap[r]))
        Rf_error("pair_wtap must contain finite values");
    }
  }

  int n_batch = as_int_scalar(n_batch_r, "n_batch");
  int batch_length = as_int_scalar(batch_length_r, "batch_length");
  int report = as_int_scalar(report_r, "report");
  int n_threads = as_int_scalar(n_threads_r, "n_threads");
  bool verbose = Rf_asLogical(verbose_r) == TRUE;

  double accept_rate = as_real_scalar(accept_rate_r, "accept_rate");
  double c0 = as_real_scalar(c0_r, "c0");
  double c1 = as_real_scalar(c1_r, "c1");
  double sigma_shape = as_real_scalar(sigma_shape_r, "sigma_shape");
  double sigma_scale = as_real_scalar(sigma_scale_r, "sigma_scale");
  double phi_lower = as_real_scalar(phi_lower_r, "phi_lower");
  double phi_upper = as_real_scalar(phi_upper_r, "phi_upper");

  if(n_batch < 1 || batch_length < 1)
    Rf_error("n_batch and batch_length must be positive");
  if(report < 1)
    report = n_batch + 1;
  if(n_threads < 1)
    n_threads = 1;

  const double *y_B = REAL(y_B_r);
  const int *trials = INTEGER(trials_r);
  const double *offset = REAL(offset_r);
  for(int i = 0; i < n_b; i++){
    if(!R_finite(offset[i]))
      Rf_error("offset must contain finite values");
    if(!R_finite(y_B[i]) || y_B[i] < 0.0 || std::floor(y_B[i]) != y_B[i])
      Rf_error("Polya-Gamma responses must be nonnegative integers");
    if(family_code == 1){
      if(trials[i] < 1)
        Rf_error("trials must contain positive integers");
      if(y_B[i] > static_cast<double>(trials[i]))
        Rf_error("Binomial responses cannot exceed trials");
    }
  }

  int theta_d = spatial ? 2 : 0;
  int n_save = n_batch * batch_length;

  SEXP theta_samples_r, beta_samples_r, omega_samples_r, eta_samples_r;
  SEXP lp_samples_r, accept_r, batch_accept_rate_r, sigma_sq_m_r;
  SEXP Sigma0_r, proposal_cov_r, out_r, names_r;

  PROTECT(theta_samples_r = Rf_allocMatrix(REALSXP, n_save, spatial ? 2 : 0));
  PROTECT(beta_samples_r = Rf_allocMatrix(REALSXP, n_save, p));
  PROTECT(omega_samples_r = Rf_allocMatrix(REALSXP, n_save, n_b));
  PROTECT(eta_samples_r = Rf_allocMatrix(REALSXP, n_save, n_b));
  PROTECT(lp_samples_r = Rf_allocVector(REALSXP, n_save));
  PROTECT(accept_r = Rf_allocMatrix(INTSXP, n_batch, spatial ? batch_length : 0));
  PROTECT(batch_accept_rate_r = Rf_allocVector(REALSXP, n_batch));
  PROTECT(sigma_sq_m_r = Rf_allocVector(REALSXP, 1));
  PROTECT(Sigma0_r = Rf_allocMatrix(REALSXP, theta_d, theta_d));
  PROTECT(proposal_cov_r = Rf_allocMatrix(REALSXP, theta_d, theta_d));

  double *theta_samples = REAL(theta_samples_r);
  double *beta_samples = REAL(beta_samples_r);
  double *omega_samples = REAL(omega_samples_r);
  double *eta_samples = REAL(eta_samples_r);
  double *lp_samples = REAL(lp_samples_r);
  int *accept = INTEGER(accept_r);
  double *batch_accept_rate = REAL(batch_accept_rate_r);
  double *Sigma0_out = REAL(Sigma0_r);
  double *proposal_cov = REAL(proposal_cov_r);

  for(int i = 0; i < n_batch * (spatial ? batch_length : 0); i++)
    accept[i] = 0;

  std::vector<double> V_beta_inv(REAL(V_beta_r), REAL(V_beta_r) + p * p);
  if(!chol_inverse_upper_in_place(V_beta_inv, p))
    Rf_error("V_beta must be positive definite");

  std::vector<double> beta(p, 0.0);
  std::vector<double> omega_B(n_b, 0.0);
  std::vector<double> alpha_mean(p + (spatial ? n_b : 0), 0.0);
  std::vector<double> pg_w(n_b, 1.0);
  std::vector<double> z(n_b, 0.0);
  std::vector<double> C_B(static_cast<R_xlen_t>(n_b) * n_b, 0.0);
  std::vector<double> C_B_inv(static_cast<R_xlen_t>(n_b) * n_b, 0.0);
  std::vector<double> theta(theta_d, 0.0);
  std::vector<double> theta_prop(theta_d, 0.0);
  std::vector<double> Sigma0(theta_d * theta_d, 0.0);
  std::vector<double> Sigma_hat(theta_d * theta_d, 0.0);
  std::vector<double> R_prop(theta_d * theta_d, 0.0);
  std::vector<double> theta_batch(batch_length * (theta_d > 0 ? theta_d : 1), 0.0);
  std::vector<double> theta_z(theta_d, 0.0);
  std::vector<double> log_vec(n_b, 0.0);

  for(int j = 0; j < p; j++)
    beta[j] = REAL(starting_beta_r)[j];

  double sigma_sq = 1.0;
  double phi = 1.0;
  double lp = 0.0;
  double sigma_sq_m = (theta_d > 0) ? 2.4 * 2.4 / static_cast<double>(theta_d) : NA_REAL;

  if(spatial){
    sigma_sq = as_real_scalar(starting_sigma_sq_r, "starting_sigma_sq");
    phi = as_real_scalar(starting_phi_r, "starting_phi");
    if(!R_finite(sigma_sq) || sigma_sq <= 0.0)
      Rf_error("starting_sigma_sq must be positive and finite");
    if(!R_finite(phi) || phi <= phi_lower || phi >= phi_upper)
      Rf_error("starting_phi must be inside (phi_lower, phi_upper)");

    theta[0] = std::log(sigma_sq);
    theta[1] = std::log((phi - phi_lower) / (phi_upper - phi));

    if(LENGTH(tuning_r) != theta_d)
      Rf_error("spatial binomial tuning must have length two");
    for(int i = 0; i < theta_d; i++)
      Sigma0[i + theta_d * i] = REAL(tuning_r)[i] * REAL(tuning_r)[i];

    lp = log_theta_binomial_spatial(
      theta.data(), omega_B.data(), n_b, pair_l, pair_k, pair_d,
      pair_wtap, n_pair, sigma_shape, sigma_scale, phi_lower, phi_upper,
      n_threads, C_B, log_vec
    );
    if(!R_finite(lp))
      Rf_error("Starting spatial binomial theta has non-finite log target density");
  }

  BayesLogit_rpg_hybrid_t pg = BayesLogit_rpg_hybrid();

  GetRNGstate();

  int save_i = 0;
  for(int b = 0; b < n_batch; b++){
    int batch_accept = 0;

    if(spatial){
      for(int i = 0; i < theta_d * theta_d; i++)
        R_prop[i] = sigma_sq_m * Sigma0[i];

      if(!chol_upper_in_place(R_prop, theta_d)){
        PutRNGstate();
        Rf_error("Proposal covariance is not positive definite");
      }
    }

    for(int j = 0; j < batch_length; j++){
      update_pg_response(y_B, trials, offset, family_code, nb_size,
                         REAL(X_B_r), n_b, p, spatial, beta, omega_B,
                         pg, pg_w, z);

      if(spatial){
        C_B_inv.assign(static_cast<R_xlen_t>(n_b) * n_b, 0.0);
        fill_C_B_from_pairs_cpp(C_B_inv.data(), n_b, pair_l, pair_k, pair_d,
                                pair_wtap, n_pair, phi, n_threads);
        if(!chol_inverse_upper_in_place(C_B_inv, n_b)){
          PutRNGstate();
          Rf_error("Current C_B is not positive definite");
        }
      }

      sample_pg_alpha(REAL(X_B_r), n_b, p, spatial, pg_w, z,
                      REAL(beta_mu_r), V_beta_inv.data(), C_B_inv,
                      sigma_sq, beta, omega_B, alpha_mean);

      if(spatial){
        for(int k = 0; k < theta_d; k++){
          double step = 0.0;
          double zk = norm_rand();
          theta_z[k] = zk;
          for(int a = 0; a <= k; a++)
            step += R_prop[a + theta_d * k] * theta_z[a];
          theta_prop[k] = theta[k] + step;
        }

        double lp_prop = log_theta_binomial_spatial(
          theta_prop.data(), omega_B.data(), n_b, pair_l, pair_k, pair_d,
          pair_wtap, n_pair, sigma_shape, sigma_scale, phi_lower, phi_upper,
          n_threads, C_B, log_vec
        );

        if(std::log(unif_rand()) < lp_prop - lp){
          theta[0] = theta_prop[0];
          theta[1] = theta_prop[1];
          lp = lp_prop;
          accept[b + n_batch * j] = 1;
          batch_accept++;
        }

        sigma_sq = std::exp(theta[0]);
        phi = phi_from_z(theta[1], phi_lower, phi_upper);

        theta_batch[j + batch_length * 0] = theta[0];
        theta_batch[j + batch_length * 1] = theta[1];
      }

      for(int k = 0; k < p; k++)
        beta_samples[save_i + n_save * k] = beta[k];

      for(int i = 0; i < n_b; i++){
        double eta = offset[i];
        for(int k = 0; k < p; k++)
          eta += REAL(X_B_r)[i + n_b * k] * beta[k];
        if(spatial)
          eta += omega_B[i];

        omega_samples[save_i + n_save * i] = omega_B[i];
        eta_samples[save_i + n_save * i] = eta;
      }

      if(spatial){
        theta_samples[save_i + n_save * 0] = theta[0];
        theta_samples[save_i + n_save * 1] = theta[1];
        lp_samples[save_i] = lp;
      } else {
        lp_samples[save_i] = NA_REAL;
      }

      save_i++;
    }

    if(spatial){
      double rhat = static_cast<double>(batch_accept) /
        static_cast<double>(batch_length);
      batch_accept_rate[b] = rhat;

      double gamma1 = 1.0 / std::pow(static_cast<double>(b + 2), c1);
      double gamma2 = c0 * gamma1;

      sigma_sq_m = std::exp(std::log(sigma_sq_m) +
        gamma2 * (rhat - accept_rate));

      if(batch_length > 1){
        covariance_from_batch(theta_batch, batch_length, theta_d, Sigma_hat);

        bool ok = true;
        for(int i = 0; i < theta_d * theta_d; i++){
          if(!R_finite(Sigma_hat[i])){
            ok = false;
            break;
          }
        }

        if(ok){
          for(int i = 0; i < theta_d * theta_d; i++)
            Sigma0[i] += gamma1 * (Sigma_hat[i] - Sigma0[i]);
        }
      }

      if(verbose && (b == 0 || ((b + 1) % report == 0))){
        Rprintf("batch %d accept = %.3f sigma_sq_m = %.4g lp = %.6g\n",
                b + 1, rhat, sigma_sq_m, lp);
      }
    } else {
      batch_accept_rate[b] = NA_REAL;
    }

    if((b + 1) % 10 == 0)
      R_CheckUserInterrupt();
  }

  PutRNGstate();

  REAL(sigma_sq_m_r)[0] = sigma_sq_m;
  for(int i = 0; i < theta_d * theta_d; i++){
    Sigma0_out[i] = Sigma0[i];
    proposal_cov[i] = sigma_sq_m * Sigma0[i];
  }

  PROTECT(out_r = Rf_allocVector(VECSXP, 10));
  SET_VECTOR_ELT(out_r, 0, theta_samples_r);
  SET_VECTOR_ELT(out_r, 1, beta_samples_r);
  SET_VECTOR_ELT(out_r, 2, omega_samples_r);
  SET_VECTOR_ELT(out_r, 3, eta_samples_r);
  SET_VECTOR_ELT(out_r, 4, lp_samples_r);
  SET_VECTOR_ELT(out_r, 5, accept_r);
  SET_VECTOR_ELT(out_r, 6, batch_accept_rate_r);
  SET_VECTOR_ELT(out_r, 7, sigma_sq_m_r);
  SET_VECTOR_ELT(out_r, 8, Sigma0_r);
  SET_VECTOR_ELT(out_r, 9, proposal_cov_r);

  PROTECT(names_r = Rf_allocVector(STRSXP, 10));
  SET_STRING_ELT(names_r, 0, Rf_mkChar("theta_z_samples"));
  SET_STRING_ELT(names_r, 1, Rf_mkChar("beta_samples"));
  SET_STRING_ELT(names_r, 2, Rf_mkChar("omega_B_samples"));
  SET_STRING_ELT(names_r, 3, Rf_mkChar("eta_B_samples"));
  SET_STRING_ELT(names_r, 4, Rf_mkChar("p.lp.samples"));
  SET_STRING_ELT(names_r, 5, Rf_mkChar("accept"));
  SET_STRING_ELT(names_r, 6, Rf_mkChar("batch.accept.rate"));
  SET_STRING_ELT(names_r, 7, Rf_mkChar("sigma_sq_m"));
  SET_STRING_ELT(names_r, 8, Rf_mkChar("Sigma0"));
  SET_STRING_ELT(names_r, 9, Rf_mkChar("proposal.cov"));
  Rf_setAttrib(out_r, R_NamesSymbol, names_r);

  UNPROTECT(12);
  return out_r;
}

/* Direct/reference construction. */

extern "C" SEXP make_C_B_tapered(SEXP plot_start_r,
                                 SEXP h_r,
                                 SEXP x_r,
                                 SEXP y_r,
                                 SEXP phi_r,
                                 SEXP gamma_r,
                                 SEXP taper_code_r,
                                 SEXP n_threads_r)
{
  int n_b, n_nz;
  double phi = as_real_scalar(phi_r, "phi");
  double gamma = as_real_scalar(gamma_r, "gamma");
  int taper_code = as_int_scalar(taper_code_r, "taper_code");
  int n_threads = as_int_scalar(n_threads_r, "n_threads");

  if(n_threads < 1)
    n_threads = 1;

  check_common_inputs(plot_start_r, h_r, x_r, y_r, phi, gamma, taper_code,
                      &n_b, &n_nz);

  int *plot_start = INTEGER(plot_start_r);
  double *h = REAL(h_r);
  double *x = REAL(x_r);
  double *y = REAL(y_r);

  SEXP C_B_r;
  PROTECT(C_B_r = Rf_allocMatrix(REALSXP, n_b, n_b));
  double *C_B = REAL(C_B_r);

  for(int i = 0; i < n_b * n_b; i++)
    C_B[i] = 0.0;

#ifdef _OPENMP
  if(n_threads > 1)
    omp_set_num_threads(n_threads);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(n_threads > 1)
#endif
  for(int l = 0; l < n_b; l++){
    int ls = plot_start[l] - 1;
    int le = plot_start[l + 1] - 1;

    for(int k = 0; k <= l; k++){
      int ks = plot_start[k] - 1;
      int ke = plot_start[k + 1] - 1;

      double val = 0.0;

      for(int a = ls; a < le; a++){
        double ha = h[a];
        double xa = x[a];
        double ya = y[a];

        for(int b = ks; b < ke; b++){
          double dx = xa - x[b];
          double dy = ya - y[b];
          double d = std::sqrt(dx * dx + dy * dy);

          if(d < gamma){
            val += ha * h[b] * std::exp(-phi * d) *
              taper_value(d, gamma, taper_code);
          }
        }
      }

      C_B[l + n_b * k] = val;
      C_B[k + n_b * l] = val;
    }

#ifndef _OPENMP
    if(l % 10 == 0)
      R_CheckUserInterrupt();
#endif
  }

  UNPROTECT(1);
  return C_B_r;
}

/* One-time precomputation of the tapered pair pattern. */

extern "C" SEXP prep_C_B_tapered_pairs(SEXP plot_start_r,
                                       SEXP h_r,
                                       SEXP x_r,
                                       SEXP y_r,
                                       SEXP gamma_r,
                                       SEXP taper_code_r,
                                       SEXP n_threads_r)
{
  int n_b, n_nz;
  double gamma = as_real_scalar(gamma_r, "gamma");
  int taper_code = as_int_scalar(taper_code_r, "taper_code");
  int n_threads = as_int_scalar(n_threads_r, "n_threads");

  if(n_threads < 1)
    n_threads = 1;

  check_common_inputs(plot_start_r, h_r, x_r, y_r, 1.0, gamma, taper_code,
                      &n_b, &n_nz);

  int *plot_start = INTEGER(plot_start_r);
  double *h = REAL(h_r);
  double *x = REAL(x_r);
  double *y = REAL(y_r);

  std::vector<int> pair_l;
  std::vector<int> pair_k;
  std::vector<double> pair_d;
  std::vector<double> pair_wtap;

#ifdef _OPENMP
  if(n_threads > 1)
    omp_set_num_threads(n_threads);
#endif

#ifdef _OPENMP
#pragma omp parallel if(n_threads > 1)
  {
    std::vector<int> local_l;
    std::vector<int> local_k;
    std::vector<double> local_d;
    std::vector<double> local_wtap;

#pragma omp for schedule(dynamic)
    for(int l = 0; l < n_b; l++){
      int ls = plot_start[l] - 1;
      int le = plot_start[l + 1] - 1;

      for(int k = 0; k <= l; k++){
        int ks = plot_start[k] - 1;
        int ke = plot_start[k + 1] - 1;

        for(int a = ls; a < le; a++){
          double ha = h[a];
          double xa = x[a];
          double ya = y[a];

          for(int b = ks; b < ke; b++){
            double dx = xa - x[b];
            double dy = ya - y[b];
            double d = std::sqrt(dx * dx + dy * dy);

            if(d < gamma){
              double tap = taper_value(d, gamma, taper_code);

              local_l.push_back(l + 1);
              local_k.push_back(k + 1);
              local_d.push_back(d);
              local_wtap.push_back(ha * h[b] * tap);
            }
          }
        }
      }
    }

#pragma omp critical
    {
      pair_l.insert(pair_l.end(), local_l.begin(), local_l.end());
      pair_k.insert(pair_k.end(), local_k.begin(), local_k.end());
      pair_d.insert(pair_d.end(), local_d.begin(), local_d.end());
      pair_wtap.insert(pair_wtap.end(), local_wtap.begin(), local_wtap.end());
    }
  }
#else
  for(int l = 0; l < n_b; l++){
    int ls = plot_start[l] - 1;
    int le = plot_start[l + 1] - 1;

    for(int k = 0; k <= l; k++){
      int ks = plot_start[k] - 1;
      int ke = plot_start[k + 1] - 1;

      for(int a = ls; a < le; a++){
        double ha = h[a];
        double xa = x[a];
        double ya = y[a];

        for(int b = ks; b < ke; b++){
          double dx = xa - x[b];
          double dy = ya - y[b];
          double d = std::sqrt(dx * dx + dy * dy);

          if(d < gamma){
            double tap = taper_value(d, gamma, taper_code);

            pair_l.push_back(l + 1);
            pair_k.push_back(k + 1);
            pair_d.push_back(d);
            pair_wtap.push_back(ha * h[b] * tap);
          }
        }
      }
    }

    if(l % 10 == 0)
      R_CheckUserInterrupt();
  }
#endif

  R_xlen_t n_pair = static_cast<R_xlen_t>(pair_l.size());

  SEXP out_r, names_r, pair_l_r, pair_k_r, pair_d_r, pair_wtap_r;
  SEXP n_b_r, gamma_out_r, taper_code_out_r;

  PROTECT(pair_l_r = Rf_allocVector(INTSXP, n_pair));
  PROTECT(pair_k_r = Rf_allocVector(INTSXP, n_pair));
  PROTECT(pair_d_r = Rf_allocVector(REALSXP, n_pair));
  PROTECT(pair_wtap_r = Rf_allocVector(REALSXP, n_pair));

  for(R_xlen_t i = 0; i < n_pair; i++){
    INTEGER(pair_l_r)[i] = pair_l[i];
    INTEGER(pair_k_r)[i] = pair_k[i];
    REAL(pair_d_r)[i] = pair_d[i];
    REAL(pair_wtap_r)[i] = pair_wtap[i];
  }

  PROTECT(n_b_r = Rf_allocVector(INTSXP, 1));
  INTEGER(n_b_r)[0] = n_b;

  PROTECT(gamma_out_r = Rf_allocVector(REALSXP, 1));
  REAL(gamma_out_r)[0] = gamma;

  PROTECT(taper_code_out_r = Rf_allocVector(INTSXP, 1));
  INTEGER(taper_code_out_r)[0] = taper_code;

  PROTECT(out_r = Rf_allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out_r, 0, pair_l_r);
  SET_VECTOR_ELT(out_r, 1, pair_k_r);
  SET_VECTOR_ELT(out_r, 2, pair_d_r);
  SET_VECTOR_ELT(out_r, 3, pair_wtap_r);
  SET_VECTOR_ELT(out_r, 4, n_b_r);
  SET_VECTOR_ELT(out_r, 5, gamma_out_r);
  SET_VECTOR_ELT(out_r, 6, taper_code_out_r);

  PROTECT(names_r = Rf_allocVector(STRSXP, 7));
  SET_STRING_ELT(names_r, 0, Rf_mkChar("pair_l"));
  SET_STRING_ELT(names_r, 1, Rf_mkChar("pair_k"));
  SET_STRING_ELT(names_r, 2, Rf_mkChar("pair_d"));
  SET_STRING_ELT(names_r, 3, Rf_mkChar("pair_wtap"));
  SET_STRING_ELT(names_r, 4, Rf_mkChar("n_b"));
  SET_STRING_ELT(names_r, 5, Rf_mkChar("gamma"));
  SET_STRING_ELT(names_r, 6, Rf_mkChar("taper_code"));

  Rf_setAttrib(out_r, R_NamesSymbol, names_r);

  UNPROTECT(9);
  return out_r;
}

/* Fast construction from precomputed pairs. */

extern "C" SEXP make_C_B_tapered_from_pairs(SEXP pair_l_r,
                                            SEXP pair_k_r,
                                            SEXP pair_d_r,
                                            SEXP pair_wtap_r,
                                            SEXP n_b_r,
                                            SEXP phi_r,
                                            SEXP n_threads_r)
{
  require_integer_vector(pair_l_r, "pair_l");
  require_integer_vector(pair_k_r, "pair_k");
  require_real_vector(pair_d_r, "pair_d");
  require_real_vector(pair_wtap_r, "pair_wtap");

  R_xlen_t n_pair = XLENGTH(pair_l_r);

  if(XLENGTH(pair_k_r) != n_pair ||
     XLENGTH(pair_d_r) != n_pair ||
     XLENGTH(pair_wtap_r) != n_pair)
    Rf_error("pair_l, pair_k, pair_d, and pair_wtap must have the same length");

  int n_b = as_int_scalar(n_b_r, "n_b");
  double phi = as_real_scalar(phi_r, "phi");
  int n_threads = as_int_scalar(n_threads_r, "n_threads");

  if(n_b <= 0)
    Rf_error("n_b must be positive");
  if(!R_finite(phi) || phi <= 0.0)
    Rf_error("phi must be finite and positive");
  if(n_threads < 1)
    n_threads = 1;

  int *pair_l = INTEGER(pair_l_r);
  int *pair_k = INTEGER(pair_k_r);
  double *pair_d = REAL(pair_d_r);
  double *pair_wtap = REAL(pair_wtap_r);

  SEXP C_B_r;
  PROTECT(C_B_r = Rf_allocMatrix(REALSXP, n_b, n_b));
  double *C_B = REAL(C_B_r);

  for(int i = 0; i < n_b * n_b; i++)
    C_B[i] = 0.0;

#ifdef _OPENMP
  if(n_threads > 1)
    omp_set_num_threads(n_threads);

#pragma omp parallel if(n_threads > 1)
  {
    std::vector<double> C_local(n_b * n_b, 0.0);

#pragma omp for schedule(static)
    for(R_xlen_t r = 0; r < n_pair; r++){
      int l = pair_l[r] - 1;
      int k = pair_k[r] - 1;

      if(l < 0 || l >= n_b || k < 0 || k >= n_b)
        continue;

      double val = pair_wtap[r] * std::exp(-phi * pair_d[r]);

      C_local[l + n_b * k] += val;

      if(l != k)
        C_local[k + n_b * l] += val;
    }

#pragma omp critical
    {
      for(int i = 0; i < n_b * n_b; i++)
        C_B[i] += C_local[i];
    }
  }
#else
  for(R_xlen_t r = 0; r < n_pair; r++){
    int l = pair_l[r] - 1;
    int k = pair_k[r] - 1;

    if(l < 0 || l >= n_b || k < 0 || k >= n_b)
      Rf_error("pair_l or pair_k contains an out-of-range plot index");

    double val = pair_wtap[r] * std::exp(-phi * pair_d[r]);

    C_B[l + n_b * k] += val;

    if(l != k)
      C_B[k + n_b * l] += val;

    if(r % 1000000 == 0)
      R_CheckUserInterrupt();
  }
#endif

  UNPROTECT(1);
  return C_B_r;
}

/* Convenience helper to build plot_start from sorted triplets. */

extern "C" SEXP make_plot_compressed_from_triplets(SEXP plot_id_r,
                                                   SEXP h_r,
                                                   SEXP x_r,
                                                   SEXP y_r,
                                                   SEXP n_b_r)
{
  require_integer_vector(plot_id_r, "plot_id");
  require_real_vector(h_r, "h");
  require_real_vector(x_r, "x");
  require_real_vector(y_r, "y");

  int n_b = as_int_scalar(n_b_r, "n_b");
  int n_nz = LENGTH(h_r);

  if(n_b <= 0)
    Rf_error("n_b must be positive");
  if(LENGTH(plot_id_r) != n_nz || LENGTH(x_r) != n_nz || LENGTH(y_r) != n_nz)
    Rf_error("plot_id, h, x, and y must have the same length");

  int *plot_id = INTEGER(plot_id_r);

  SEXP plot_start_r;
  PROTECT(plot_start_r = Rf_allocVector(INTSXP, n_b + 1));
  int *plot_start = INTEGER(plot_start_r);

  int pos = 0;
  for(int l = 1; l <= n_b; l++){
    plot_start[l - 1] = pos + 1;

    while(pos < n_nz && plot_id[pos] == l)
      pos++;

    if(pos < n_nz && plot_id[pos] < l)
      Rf_error("plot_id must be sorted in nondecreasing order");
  }

  plot_start[n_b] = n_nz + 1;

  SEXP out_r, names_r;
  PROTECT(out_r = Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out_r, 0, plot_start_r);
  SET_VECTOR_ELT(out_r, 1, h_r);
  SET_VECTOR_ELT(out_r, 2, x_r);
  SET_VECTOR_ELT(out_r, 3, y_r);

  PROTECT(names_r = Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(names_r, 0, Rf_mkChar("plot_start"));
  SET_STRING_ELT(names_r, 1, Rf_mkChar("h"));
  SET_STRING_ELT(names_r, 2, Rf_mkChar("x"));
  SET_STRING_ELT(names_r, 3, Rf_mkChar("y"));

  Rf_setAttrib(out_r, R_NamesSymbol, names_r);

  UNPROTECT(3);
  return out_r;
}



/* Fine-pixel to observed-support tapered cross-covariance.

   Returns C_AB with dimension n_a x n_b, where

     C_AB(i,l) = sum_j h_lj exp(-phi d(i,j)) T_gamma(d(i,j)).

   This is the tapered fine-to-observed-support cross-covariance used in
   conditional fine-support prediction. It is intentionally written as
   simple flat loops; the R code handles the conditional Gaussian algebra.
*/

extern "C" SEXP make_C_AB_tapered(SEXP a_x_r,
                                  SEXP a_y_r,
                                  SEXP plot_start_r,
                                  SEXP h_r,
                                  SEXP x_r,
                                  SEXP y_r,
                                  SEXP phi_r,
                                  SEXP gamma_r,
                                  SEXP taper_code_r,
                                  SEXP n_threads_r)
{
  require_real_vector(a_x_r, "a_x");
  require_real_vector(a_y_r, "a_y");

  if(XLENGTH(a_x_r) != XLENGTH(a_y_r))
    Rf_error("a_x and a_y must have the same length");

  int n_a = LENGTH(a_x_r);
  int n_b, n_nz;
  double phi = as_real_scalar(phi_r, "phi");
  double gamma = as_real_scalar(gamma_r, "gamma");
  int taper_code = as_int_scalar(taper_code_r, "taper_code");
  int n_threads = as_int_scalar(n_threads_r, "n_threads");

  if(n_threads < 1)
    n_threads = 1;

  check_common_inputs(plot_start_r, h_r, x_r, y_r, phi, gamma, taper_code,
                      &n_b, &n_nz);

  double *a_x = REAL(a_x_r);
  double *a_y = REAL(a_y_r);
  int *plot_start = INTEGER(plot_start_r);
  double *h = REAL(h_r);
  double *x = REAL(x_r);
  double *y = REAL(y_r);

  SEXP C_AB_r;
  PROTECT(C_AB_r = Rf_allocMatrix(REALSXP, n_a, n_b));
  double *C_AB = REAL(C_AB_r);

  for(int i = 0; i < n_a * n_b; i++)
    C_AB[i] = 0.0;

#ifdef _OPENMP
  if(n_threads > 1)
    omp_set_num_threads(n_threads);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(n_threads > 1)
#endif
  for(int i = 0; i < n_a; i++){
    double xi = a_x[i];
    double yi = a_y[i];

    for(int l = 0; l < n_b; l++){
      int ls = plot_start[l] - 1;
      int le = plot_start[l + 1] - 1;
      double val = 0.0;

      for(int j = ls; j < le; j++){
        double dx = xi - x[j];
        double dy = yi - y[j];
        double d = std::sqrt(dx * dx + dy * dy);

        if(d < gamma)
          val += h[j] * std::exp(-phi * d) * taper_value(d, gamma, taper_code);
      }

      C_AB[i + n_a * l] = val;
    }

#ifndef _OPENMP
    if(i % 100 == 0)
      R_CheckUserInterrupt();
#endif
  }

  UNPROTECT(1);
  return C_AB_r;
}

/* Dense tapered fine-support covariance for small development examples.

   Returns C_A with dimension n_a x n_a.  This should not be used for large
   grids. It is provided only for optional sampling of omega_A from the
   full conditional covariance. The default fine-support prediction workflow
   uses the conditional mean and does not need C_A.
*/

extern "C" SEXP make_C_A_tapered_dense(SEXP a_x_r,
                                       SEXP a_y_r,
                                       SEXP phi_r,
                                       SEXP gamma_r,
                                       SEXP taper_code_r,
                                       SEXP n_threads_r)
{
  require_real_vector(a_x_r, "a_x");
  require_real_vector(a_y_r, "a_y");

  if(XLENGTH(a_x_r) != XLENGTH(a_y_r))
    Rf_error("a_x and a_y must have the same length");

  int n_a = LENGTH(a_x_r);
  double phi = as_real_scalar(phi_r, "phi");
  double gamma = as_real_scalar(gamma_r, "gamma");
  int taper_code = as_int_scalar(taper_code_r, "taper_code");
  int n_threads = as_int_scalar(n_threads_r, "n_threads");

  if(!R_finite(phi) || phi <= 0.0)
    Rf_error("phi must be finite and positive");
  if(!R_finite(gamma) || gamma <= 0.0)
    Rf_error("gamma must be finite and positive");
  if(taper_code != 1 && taper_code != 2)
    Rf_error("taper_code must be 1 (Wendland) or 2 (spherical)");
  if(n_threads < 1)
    n_threads = 1;

  R_xlen_t n_elem = static_cast<R_xlen_t>(n_a) * static_cast<R_xlen_t>(n_a);
  if(n_elem > static_cast<R_xlen_t>(INT_MAX))
    Rf_error("Dense fine-support covariance is too large for this implementation; use fewer prediction cells or method = 'mean'");

  double *a_x = REAL(a_x_r);
  double *a_y = REAL(a_y_r);

  SEXP C_A_r;
  PROTECT(C_A_r = Rf_allocMatrix(REALSXP, n_a, n_a));
  double *C_A = REAL(C_A_r);

  for(R_xlen_t i = 0; i < n_elem; i++)
    C_A[i] = 0.0;

#ifdef _OPENMP
  if(n_threads > 1)
    omp_set_num_threads(n_threads);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(n_threads > 1)
#endif
  for(int i = 0; i < n_a; i++){
    C_A[static_cast<R_xlen_t>(i) + static_cast<R_xlen_t>(n_a) * i] = 1.0;

    for(int j = 0; j < i; j++){
      double dx = a_x[i] - a_x[j];
      double dy = a_y[i] - a_y[j];
      double d = std::sqrt(dx * dx + dy * dy);
      double val = 0.0;

      if(d < gamma)
        val = std::exp(-phi * d) * taper_value(d, gamma, taper_code);

      C_A[static_cast<R_xlen_t>(i) + static_cast<R_xlen_t>(n_a) * j] = val;
      C_A[static_cast<R_xlen_t>(j) + static_cast<R_xlen_t>(n_a) * i] = val;
    }

#ifndef _OPENMP
    if(i % 100 == 0)
      R_CheckUserInterrupt();
#endif
  }

  UNPROTECT(1);
  return C_A_r;
}

static const R_CallMethodDef CallEntries[] = {
  {"gaussian_metrop_sampler", (DL_FUNC) &gaussian_metrop_sampler, 27},
  {"binomial_pg_sampler", (DL_FUNC) &binomial_pg_sampler, 30},
  {"make_C_B_tapered", (DL_FUNC) &make_C_B_tapered, 8},
  {"prep_C_B_tapered_pairs", (DL_FUNC) &prep_C_B_tapered_pairs, 7},
  {"make_C_B_tapered_from_pairs", (DL_FUNC) &make_C_B_tapered_from_pairs, 7},
  {"make_C_AB_tapered", (DL_FUNC) &make_C_AB_tapered, 10},
  {"make_C_A_tapered_dense", (DL_FUNC) &make_C_A_tapered_dense, 6},
  {"make_plot_compressed_from_triplets", (DL_FUNC) &make_plot_compressed_from_triplets, 5},
  {NULL, NULL, 0}
};

extern "C" void R_init_cosTapered(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
