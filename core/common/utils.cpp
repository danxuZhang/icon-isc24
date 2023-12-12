#include "utils.hpp"
#include <iostream>

void utils_muphys::calc_dz(array_1d_t<real_t> &z, array_1d_t<real_t> &dz,
                           size_t &ncells, size_t &nlev) {
  dz.resize(ncells * nlev);
  array_2d_t<real_t> zh(nlev + 1, array_1d_t<real_t>(ncells));

  for (size_t i = 0; i < ncells; i++) {
    zh[nlev][i] = (static_cast<real_t>(3.0) * z[i + (nlev - 1) * (ncells)] -
                   z[i + (nlev - 2) * (ncells)]) *
                  static_cast<real_t>(0.5);
  }

  // The loop is intentionally i<nlev; since we are using an unsigned integer
  // data type, when i reaches 0, and you try to decrement further, (to -1), it
  // wraps to the maximum value representable by size_t.
  for (size_t i = nlev - 1; i < nlev; --i) {
    for (size_t j = 0; j < ncells; j++) {
      zh[i][j] = static_cast<real_t>(2.0) * z[j + (i * ncells)] - zh[i + 1][j];
      dz[i * ncells + j] = -zh[i + 1][j] + zh[i][j];
    }
  }
}