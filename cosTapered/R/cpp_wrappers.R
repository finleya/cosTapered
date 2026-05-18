make_C_B_tapered_Rcall <- function(plot_start, h, x, y, phi, gamma,
                                   taper_code = 1L, n_threads = 1L) {
  .Call("make_C_B_tapered",
        as.integer(plot_start),
        as.numeric(h),
        as.numeric(x),
        as.numeric(y),
        as.numeric(phi),
        as.numeric(gamma),
        as.integer(taper_code),
        as.integer(n_threads),
        PACKAGE = "cosTapered")
}

prep_C_B_tapered_pairs_Rcall <- function(plot_start, h, x, y, gamma,
                                          taper_code = 1L, n_threads = 1L) {
  .Call("prep_C_B_tapered_pairs",
        as.integer(plot_start),
        as.numeric(h),
        as.numeric(x),
        as.numeric(y),
        as.numeric(gamma),
        as.integer(taper_code),
        as.integer(n_threads),
        PACKAGE = "cosTapered")
}

make_C_B_tapered_from_pairs_Rcall <- function(C_B_pairs, phi,
                                              n_threads = 1L) {
  .Call("make_C_B_tapered_from_pairs",
        as.integer(C_B_pairs$pair_l),
        as.integer(C_B_pairs$pair_k),
        as.numeric(C_B_pairs$pair_d),
        as.numeric(C_B_pairs$pair_wtap),
        as.integer(C_B_pairs$n_b),
        as.numeric(phi),
        as.integer(n_threads),
        PACKAGE = "cosTapered")
}

make_C_pred_B_tapered_Rcall <- function(pred_coords, H_BA_comp, phi, gamma,
                                        taper_code = 1L, n_threads = 1L) {
  pred_coords <- as.matrix(pred_coords)

  .Call("make_C_AB_tapered",
        as.numeric(pred_coords[, 1]),
        as.numeric(pred_coords[, 2]),
        as.integer(H_BA_comp$plot_start),
        as.numeric(H_BA_comp$h),
        as.numeric(H_BA_comp$x),
        as.numeric(H_BA_comp$y),
        as.numeric(phi),
        as.numeric(gamma),
        as.integer(taper_code),
        as.integer(n_threads),
        PACKAGE = "cosTapered")
}

make_C_pred_tapered_dense_Rcall <- function(pred_coords, phi, gamma,
                                            taper_code = 1L, n_threads = 1L) {
  pred_coords <- as.matrix(pred_coords)

  .Call("make_C_A_tapered_dense",
        as.numeric(pred_coords[, 1]),
        as.numeric(pred_coords[, 2]),
        as.numeric(phi),
        as.numeric(gamma),
        as.integer(taper_code),
        as.integer(n_threads),
        PACKAGE = "cosTapered")
}
