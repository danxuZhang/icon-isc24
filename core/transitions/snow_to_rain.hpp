#pragma once

#include "../common/constants.hpp"
#include "../common/types.hpp"
#include <cmath>

namespace transition {

constexpr real_t c1_sr = 79.6863;     // Constants in melting formula
constexpr real_t c2_sr = 0.612654E-3; // Constants in melting formula
constexpr real_t a_sr =
    graupel_ct::tx - static_cast<real_t>(389.5); // melting prefactor
constexpr real_t b_sr =
    static_cast<real_t>(4.0) / static_cast<real_t>(5.0); // melting exponent

/**
 * @brief Melting of snow to form rain
 *
 * @param [in] t Temperature
 * @param [in] p Ambient pressure
 * @param [in] rho Ambient density
 * @param [in] dvsw0  qv-qsat_water(T0)
 * @param [in] qs Snow specific mass
 * @return conversion rate from snow to rain
 */
TARGET real_t snow_to_rain(real_t t, real_t p, real_t rho, real_t dvsw0,
                           real_t qs) {

  return (t > fmax(thermodyn::tmelt,
                   thermodyn::tmelt - graupel_ct::tx * dvsw0) &&
          qs > graupel_ct::qmin)
             ? (c1_sr / p + c2_sr) * (t - thermodyn::tmelt + a_sr * dvsw0) *
                   pow(qs * rho, b_sr)
             : static_cast<real_t>(0.0);
}
} // namespace transition
