#include <R.h>
#include <Rinternals.h>
#include <cmath>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

// Build the self-covariance cache using each unordered pair once. Sorting
// by y lets the inner loop stop as soon as the vertical distance exceeds
// the taper. Per-thread accumulators avoid locks in the pair loop.
struct GridCell {
  int x;
  int y;
  double h;
};

extern "C" SEXP cos_areal_grid_self_weights_cpp(SEXP x_r, SEXP y_r, SEXP h_r,
                                           SEXP n_bin_r, SEXP n_threads_r)
{
  if (TYPEOF(x_r) != INTSXP || TYPEOF(y_r) != INTSXP || TYPEOF(h_r) != REALSXP)
    Rf_error("Invalid lattice covariance inputs");
  const int n = Rf_length(h_r), nbin = Rf_asInteger(n_bin_r);
  if (n < 1 || nbin < 1 || nbin > 1000000 ||
      Rf_length(x_r) != n || Rf_length(y_r) != n)
    Rf_error("Incompatible lattice covariance dimensions");
  // Limit temporary accumulators to 64 MiB even if many threads are requested.
  int nt = std::min(n, std::max(1, Rf_asInteger(n_threads_r)));
  nt = std::min(nt, (int)((64 * 1024 * 1024) / (sizeof(double) * nbin)));
#ifndef _OPENMP
  nt = 1;
#endif
  const int *x = INTEGER(x_r), *y = INTEGER(y_r);
  const double *h = REAL(h_r);
  GridCell *cells = (GridCell *)R_alloc(n, sizeof(GridCell));
  for (int i = 0; i < n; ++i) {
    if (x[i] == NA_INTEGER || y[i] == NA_INTEGER || !R_finite(h[i]) || h[i] < 0)
      Rf_error("Cell coordinates and nonnegative weights must be finite");
    cells[i] = {x[i], y[i], h[i]};
  }
  std::sort(cells, cells + n, [](const GridCell &a, const GridCell &b) {
    return a.y < b.y || (a.y == b.y && a.x < b.x);
  });
  const int radius = (int)std::ceil(std::sqrt((double)nbin));
  double *scratch = (double *)R_alloc((size_t)nt * nbin, sizeof(double));
  std::fill(scratch, scratch + (size_t)nt * nbin, 0.0);
  // R-managed temporary storage is reclaimed if a user interrupts a batch.
  for (int first = 0; first < n; ) {
    const int end = first + std::min(4096, n - first);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
#endif
    for (int i = first; i < end; ++i) {
      int worker = 0;
#ifdef _OPENMP
      worker = omp_get_thread_num();
#endif
      double *weights = scratch + (size_t)worker * nbin;
      const GridCell &a = cells[i];
      weights[0] += a.h * a.h;
      for (int j = i + 1; j < n; ++j) {
        const long long dy = (long long)cells[j].y - a.y;
        if (dy >= radius) break;
        const long long dx = (long long)cells[j].x - a.x;
        if (dx <= -radius || dx >= radius) continue;
        const long long d2 = dx * dx + dy * dy;
        // Both (i,j) and (j,i) contribute to the original double sum.
        if (d2 < nbin) weights[d2] += 2.0 * a.h * cells[j].h;
      }
    }
    first = end;
    R_CheckUserInterrupt();
  }
  SEXP out = PROTECT(Rf_allocVector(REALSXP, nbin));
  std::fill(REAL(out), REAL(out) + nbin, 0.0);
  for (int worker = 0; worker < nt; ++worker)
    for (int k = 0; k < nbin; ++k)
      REAL(out)[k] += scratch[(size_t)worker * nbin + k];
  UNPROTECT(1);
  return out;
}

// Sum support-pair weights by exact squared lattice distance. All geometry
// work is done once, independent of posterior phi. Columns are observed
// supports, so even a single forest target can use all OpenMP workers.
extern "C" SEXP cos_areal_grid_cross_weights(SEXP ux_r, SEXP uy_r, SEXP uh_r,
                                            SEXP bx_r, SEXP by_r, SEXP bh_r,
                                            SEXP start_r, SEXP n_bin_r,
                                            SEXP n_threads_r)
{
  if (TYPEOF(ux_r) != INTSXP || TYPEOF(uy_r) != INTSXP ||
      TYPEOF(bx_r) != INTSXP || TYPEOF(by_r) != INTSXP ||
      TYPEOF(start_r) != INTSXP || TYPEOF(uh_r) != REALSXP ||
      TYPEOF(bh_r) != REALSXP)
    Rf_error("Invalid lattice covariance inputs");
  const int nu = Rf_length(uh_r), nbcell = Rf_length(bh_r);
  const int nb = Rf_length(start_r) - 1;
  const int nbin = Rf_asInteger(n_bin_r);
  const int nt = std::max(1, Rf_asInteger(n_threads_r));
  if (nb < 1 || nbin < 1 || nbin > 1000000 ||
      Rf_length(ux_r) != nu || Rf_length(uy_r) != nu ||
      Rf_length(bx_r) != nbcell || Rf_length(by_r) != nbcell)
    Rf_error("Incompatible lattice covariance dimensions");
  const int *ux = INTEGER(ux_r), *uy = INTEGER(uy_r);
  const int *bx = INTEGER(bx_r), *by = INTEGER(by_r), *start = INTEGER(start_r);
  const double *uh = REAL(uh_r), *bh = REAL(bh_r);
  if (start[0] != 1 || start[nb] != nbcell + 1)
    Rf_error("Invalid observed-support offsets");
  for (int k = 0; k < nb; ++k)
    if (start[k] < 1 || start[k + 1] < start[k])
      Rf_error("Invalid observed-support offsets");
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, nbin, nb));
  double *weights = REAL(out);
  std::fill(weights, weights + (R_xlen_t)nbin * nb, 0.0);
  const int radius = (int)std::ceil(std::sqrt((double)nbin));
  // Batches bound the time between user-interrupt checks. No R API calls
  // occur inside the parallel region. Each thread owns an output column.
  for (int first = 0; first < nb; first += nt) {
    const int end = std::min(nb, first + nt);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if(nt > 1)
#endif
    for (int k = first; k < end; ++k) {
      double *wk = weights + (R_xlen_t)nbin * k;
      for (int j = start[k] - 1; j < start[k + 1] - 1; ++j) {
        for (int i = 0; i < nu; ++i) {
          const long long dx = (long long)ux[i] - bx[j];
          const long long dy = (long long)uy[i] - by[j];
          if (dx <= -radius || dx >= radius || dy <= -radius || dy >= radius)
            continue;
          const long long d2 = dx * dx + dy * dy;
          if (d2 < nbin) wk[d2] += uh[i] * bh[j];
        }
      }
    }
    R_CheckUserInterrupt();
  }
  UNPROTECT(1);
  return out;
}
