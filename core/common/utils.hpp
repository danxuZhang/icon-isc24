#pragma once
#include "types.hpp"

namespace utils_muphys {
void calc_dz(array_1d_t<real_t> &z, array_1d_t<real_t> &dz, size_t &ncells,
             size_t &nlev);
}
