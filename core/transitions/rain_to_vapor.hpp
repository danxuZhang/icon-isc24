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

constexpr real_t b1_rv =
    0.16667; // exponent in power-law relation for mass density
constexpr real_t b2_rv = 0.55555;
constexpr real_t c1_rv = 0.61;
constexpr real_t c2_rv = -0.0163;
constexpr real_t c3_rv = 1.111e-4;
constexpr real_t a1_rv = 1.536e-3;
constexpr real_t a2_rv = 1.0E+0; // constant in rain evap formula
constexpr real_t a3_rv =
    19.0621E+0; // prefactor (from gamma dist. and properties of air/water)

/**
* @brief Calculates the mass transfer from rain to water vapor (evaporation)
*
* @param [in] t Temperature (K)
* @param [in] rho Ambient density (kg/m^3)
* @param [in] qc Specific humidity of cloud (kg/kg)
* @param [in] qr Specific humidity of rain (kg/kg)
* @param [in] dvsw qv - qsat_water(T) (kg/kg)
* @param [in] dt Time step (s)
* @return The mass transfer rate from rain to water vapor (kg/kg/s)
*/
TARGET real_t rain_to_vapor(real_t t, real_t rho, real_t qc, real_t qr,
                            real_t dvsw, real_t dt) {

  if (qr > graupel_ct::qmin && (dvsw + qc <= 0.0)) {
    real_t tc = t - thermodyn::tmelt;
    real_t evap_max = (c1_rv + tc * (c2_rv + c3_rv * tc)) * (-dvsw) / dt;
    return fmin(a1_rv * (a2_rv + a3_rv * pow((qr * rho), b1_rv)) * (-dvsw) *
                    pow((qr * rho), b2_rv),
                evap_max);
  }

  return static_cast<real_t>(0.0);
}

} // namespace transition
