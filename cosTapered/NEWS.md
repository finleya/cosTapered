# cosTapered 0.0.3

* Accelerate areal prediction on square raster grids by caching covariance
  weights at exact squared lattice distances. Use parallel C++ pair loops,
  symmetry, and taper bounds for self-covariances, and parallel preparation of
  covariances with observed supports. Preserve Gaussian, binomial, and negative-binomial
  prediction semantics, including marginal spatial uncertainty.
* Record raster grid metadata in newly constructed areal blocks. Keep the
  direct calculation for older block objects, incompatible grids, small
  supports, and calculations exceeding internal memory limits.
