#pragma once

#include "../common/constants.hpp"
#include "../common/types.hpp"
#include <cmath>

namespace property {

constexpr real_t mi_max = 1.0e-09; // maximum mass of cloud ice crystals

/**
 * @brief TODO
 * @param [in] qi Ice specific mass
 * @param [in] ni Ice crystal number
 * @return ice mass
 */
TARGET real_t ice_mass(real_t qi, real_t ni) {
  return fmax(graupel_ct::m0_ice, fmin(qi / ni, mi_max));
}
} // namespace property
