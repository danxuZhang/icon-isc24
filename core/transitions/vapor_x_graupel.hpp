#pragma once

#include "../common/constants.hpp"
#include "../common/types.hpp"
#include <cmath>

namespace transition {

constexpr real_t a1_vg = 0.398561;
constexpr real_t a2_vg = -0.00152398;
constexpr real_t a3 = 2554.99;
constexpr real_t a4 = 2.6531E-7;
constexpr real_t a5 = 0.153907;
constexpr real_t a6 = -7.86703e-07;
constexpr real_t a7 = 0.0418521;
constexpr real_t a8 = -4.7524E-8;
constexpr real_t b_vg = 0.6;

/**
 * @brief Graupel-vapor exchange rate
 *
 * @param [in] t Temperature
 * @param [in] p Ambient pressure
 * @param [in] rho Ambient density
 * @param [in] qg Graupel specific mass
 * @param [in] dvsw qv-qsat_water(T)
 * @param [in] dvsi qv-qsat_ice(T)
 * @param [in] dvsw0 qv-qsat_water(T0)
 * @param [in] dt Time step
 * @return  TODO
 */
TARGET real_t vapor_x_graupel(real_t t, real_t p, real_t rho, real_t qg,
                              real_t dvsw, real_t dvsi, real_t dvsw0,
                              real_t dt) {

  real_t result = 0.0;

  if (qg > graupel_ct::qmin) {
    if (t < thermodyn::tmelt) {
      result =
          (a1_vg + a2_vg * t + a3 / p + a4 * p) * dvsi * pow(qg * rho, b_vg);
    } else {
      if (t > (thermodyn::tmelt - graupel_ct::tx * dvsw0)) {
        result = (a5 + a6 * p) * fmin(static_cast<real_t>(0.0), dvsw0) *
                 pow(qg * rho, b_vg);
      } else {
        result = (a7 + a8 * p) * dvsw * pow(qg * rho, b_vg);
      }
    }
    result = fmax(result, -qg / dt);
  }

  return result;
}
} // namespace transition
