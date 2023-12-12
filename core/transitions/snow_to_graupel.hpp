#pragma once

#include "../common/constants.hpp"
#include "../common/types.hpp"
#include <cmath>

namespace transition {

constexpr real_t a_rim_ct = 0.5; /// Constants in riming formula
constexpr real_t b_rim_ct = static_cast<real_t>(3.0) / static_cast<real_t>(4.0);

/**
 * @brief TODO
 * @param [in] t ambient temperature
 * @param [in] rho ambient density
 * @param [in] qc cloud specific mass
 * @param [in] qs snow specific mass
 * @returns convertion rate
 */
TARGET real_t snow_to_graupel(real_t t, real_t rho, real_t qc, real_t qs) {

  return (fmin(qc, qs) > graupel_ct::qmin && t > graupel_ct::tfrz_hom)
             ? a_rim_ct * qc * pow(qs * rho, b_rim_ct)
             : static_cast<real_t>(0.0);
}

} // namespace transition
