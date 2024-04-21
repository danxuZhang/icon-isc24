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

constexpr real_t ami = 130.0; // Formfactor for mass-size relation of cld ice
constexpr real_t b_exp =
    -0.67; // exp. for conv. (-1 + 0.33) of ice mass to sfc area

/**
 * @brief Calculates the rate of vapor deposition to ice
 *
 * This function calculates the rate at which water vapor is deposited onto ice crystals,
 * or the rate at which ice crystals sublimate, depending on the ambient conditions and
 * the state of the ice. The deposition/sublimation rate is based on the specific humidity
 * of ice, ice crystal mass, deposition factor, vapor excess with respect to ice saturation,
 * and the time step.
 *
 * @param [in] qi Specific humidity of ice (kg/kg)
 * @param [in] mi Ice crystal mass (kg)
 * @param [in] eta Deposition factor (dimensionless)
 * @param [in] dvsi Vapor excess with respect to ice saturation (kg/kg)
 * @param [in] dt Time step (s)
 * @return The rate of vapor deposition to ice (kg/kg/s)
 */

TARGET real_t vapor_x_ice(real_t qi, real_t mi, real_t eta, real_t dvsi,
                          real_t dt) {
  const real_t a_fact =
      static_cast<real_t>(4.0) *
      pow(ami, static_cast<real_t>(-1.0) / static_cast<real_t>(3.0));
  real_t result = 0.0;

  if (qi > graupel_ct::qmin) {
    result = (a_fact * eta) * qi * pow(mi, b_exp) * dvsi;

    if (result > 0.) {
      result = fmin(result, dvsi / dt);
    } else {
      result = fmax(result, dvsi / dt);
      result = fmax(result, -qi / dt);
    }
  }

  return result;
}
} // namespace transition
