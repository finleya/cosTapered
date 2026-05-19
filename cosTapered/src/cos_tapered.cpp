#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

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
