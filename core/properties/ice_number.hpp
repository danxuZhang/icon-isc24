#pragma once

#include "../common/constants.hpp"
#include "../common/types.hpp"
#include <cmath>

namespace property {

constexpr real_t a_coop = 5.000;  // parameter in cooper fit
constexpr real_t b_coop = 0.304;  // parameter in cooper fit
constexpr real_t nimax = 250.e+3; // maximal number of ice crystals

/**
 * @brief Ice number following cooper
 *
 * @param [in] t Ambient temperature (kelvin)
 * @param [in] rho Ambient density
 * @return Ice number
 */
TARGET real_t ice_number(real_t t, real_t rho) {
  return fmin(nimax, a_coop * exp(b_coop * (thermodyn::tmelt - t))) / rho;
}

} // namespace property
