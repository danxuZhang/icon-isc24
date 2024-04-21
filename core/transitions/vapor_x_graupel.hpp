// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------
//
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
* @brief Calculates the graupel-vapor exchange rate
*
* This function calculates the exchange rate between graupel and water vapor,
* which can be positive (deposition) or negative (sublimation), depending on
* the ambient conditions and the state of the graupel.
*
* @param [in] t Temperature (K)
* @param [in] p Ambient pressure (Pa)
* @param [in] rho Ambient density (kg/m^3)
* @param [in] qg Graupel specific mass (kg/kg)
* @param [in] dvsw qv - qsat_water(T) (kg/kg)
* @param [in] dvsi qv - qsat_ice(T) (kg/kg)
* @param [in] dvsw0 qv - qsat_water(T0) (kg/kg)
* @param [in] dt Time step (s)
* @return The graupel-vapor exchange rate (kg/kg/s)
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
