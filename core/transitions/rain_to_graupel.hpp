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

constexpr real_t tfrz_rain = thermodyn::tmelt - static_cast<real_t>(2.0);
constexpr real_t a1 =
    9.95e-5; // FR: 1. coefficient for immersion raindrop freezing: alpha_if
constexpr real_t b1 =
    static_cast<real_t>(7.0) /
    static_cast<real_t>(
        4.0); // FR: 2. coefficient for immersion raindrop freezing: a_if
constexpr real_t c1 = 1.68; //  coefficient for raindrop freezing
constexpr real_t c2 =
    0.66; // FR: 2. coefficient for immersion raindrop freezing: a_if
constexpr real_t c3 =
    1.0; // FR: 2. coefficient for immersion raindrop freezing: a_if
constexpr real_t c4 =
    0.1; // FR: 2. coefficient for immersion raindrop freezing: a_if
constexpr real_t a2 = 1.24E-3; //  (PI/24)*EIR*V0R*Gamma(6.5)*AR**(-5/8)
constexpr real_t b2 = static_cast<real_t>(13.0) / static_cast<real_t>(8.0);
constexpr real_t qs_crit = 1.e-7;

/**
 * @brief Freezing rain
 *
 * @param [in] t Temperature
 * @param [in] rho Ambient density
 * @param [in] qr Specific humidity of rain
 * @param [in] qc Cloud liquid specific mass
 * @param [in] qi Cloud ice specific mass
 * @param [in] qs Snow specific mass
 * @param [in] mi Ice crystal mass
 * @param [in] dvsw qv-qsat_water(T)
 * @param [in] dt Time step
 * @return convertion rate from graupel to rain
 */
TARGET real_t rain_to_graupel(real_t t, real_t rho, real_t qc, real_t qr,
                              real_t qi, real_t qs, real_t mi, real_t dvsw,
                              real_t dt) {

  real_t result = 0.0;

  if (qr > graupel_ct::qmin && t < tfrz_rain) {
    if (t > graupel_ct::tfrz_hom) {
      if (dvsw + qc <= 0.0 || qr > c4 * qc) {
        result = (exp(c2 * (tfrz_rain - t)) - c3) * (a1 * pow((qr * rho), b1));
      }
    } else {
      result = qr / dt;
    }
  }

  if (fmin(qi, qr) > graupel_ct::qmin &&
      qs > qs_crit) { // ! rain + ice creating graupel
    result = result + a2 * (qi / mi) * pow(qr, b2);
  }

  return result;
}

} // namespace transition
