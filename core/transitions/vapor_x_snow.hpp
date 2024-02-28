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

constexpr real_t nu = 1.75e-5; // kinematic viscosity of air
constexpr real_t a0_vs = 1.0;
constexpr real_t a2_vs = -(graupel_ct::v1s + 1.0) / static_cast<real_t>(2.0);
constexpr real_t eps = 1.e-15;
constexpr real_t qs_lim = 1.e-7;
constexpr real_t cnx = 4.0;
constexpr real_t b_vs = 0.8;
constexpr real_t c1_vs = 31282.3;
constexpr real_t c2_vs = 0.241897;
constexpr real_t c3_vs = 0.28003;
constexpr real_t c4_vs = -0.146293E-6;

/**
 * @brief TODO
 * @param [in] t Temperature
 * @param [in] p Ambient pressure
 * @param [in] rho Ambient density
 * @param [in] qs Snow specific mass
 * @param [in] ns Snow number
 * @param [in] lambda  Slope parameter (lambda) snow
 * @param [in] eta Deposition factor
 * @param [in] ice_dep Limiter for vapor dep on snow
 * @param [in] dvsw qv-qsat_water(T)
 * @param [in] dvsi  qv-qsat_ice(T)
 * @param [in] dvsw0 qv-qsat_water(T0)
 * @param [in] dt Time step
 * @return Rate of vapor deposition to snow
 */
TARGET real_t vapor_x_snow(real_t t, real_t p, real_t rho, real_t qs, real_t ns,
                           real_t lambda, real_t eta, real_t ice_dep,
                           real_t dvsw, real_t dvsi, real_t dvsw0, real_t dt) {
  const real_t a1_vs = static_cast<real_t>(0.4182) * sqrt(graupel_ct::v0s / nu);
  real_t result = 0.0;

  if (qs > graupel_ct::qmin) {
    if (t < thermodyn::tmelt) {
      result = (cnx * ns * eta / rho) * (a0_vs + a1_vs * pow(lambda, a2_vs)) *
               dvsi / (lambda * lambda + eps);

      // GZ: This limitation, which was missing in the original graupel scheme,
      // is crucial for numerical stability in the tropics!
      // a meaningful distiction between cloud ice and snow

      if (result > 0.0)
        result = fmin(result, dvsi / dt - ice_dep);
      if (qs <= qs_lim)
        result = fmin(result, static_cast<real_t>(0.0));
    } else {
      if (t > (thermodyn::tmelt - graupel_ct::tx * dvsw0)) {
        result = (c1_vs / p + c2_vs) * fmin(static_cast<real_t>(0.0), dvsw0) *
                 pow(qs * rho, b_vs);
      } else {
        result = (c3_vs + c4_vs * p) * dvsw * pow(qs * rho, b_vs);
      }
    }
    result = fmax(result, -qs / dt);
  }

  return result;
}

} // namespace transition
