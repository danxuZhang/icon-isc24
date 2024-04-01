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

/**
 * @brief Calculates the homogeneous freezing rate between cloud water and ice.
 *
 * This function calculates the rate at which cloud water is converted to ice through the
 * homogeneous freezing process, or the rate at which ice is converted back to cloud water
 * through melting. The freezing/melting rate is based on the temperature, cloud water
 * specific mass, ice specific mass, and the time step.
 *
 * If the cloud water specific mass is above a threshold (graupel_ct::qmin) and the temperature
 * is below the homogeneous freezing threshold (graupel_ct::tfrz_hom), the freezing rate is
 * calculated as the cloud water specific mass divided by the time step.
 *
 * If the ice specific mass is above a threshold (graupel_ct::qmin) and the temperature is
 * above the melting point (thermodyn::tmelt), the melting rate is calculated as the negative
 * of the ice specific mass divided by the time step.
 *
 * @param [in] t Temperature (K).
 * @param [in] qc Cloud water specific mass (kg/kg).
 * @param [in] qi Ice specific mass (kg/kg).
 * @param [in] dt Time step (s).
 * @return The homogeneous freezing rate (positive for freezing, negative for melting) (kg/kg/s).
 */
TARGET real_t cloud_x_ice(real_t t, real_t qc, real_t qi, real_t dt) {
  real_t result = 0.0;

  if (qc > graupel_ct::qmin && t < graupel_ct::tfrz_hom)
    result = qc / dt;

  if (qi > graupel_ct::qmin && t > thermodyn::tmelt)
    result = -qi / dt;

  return result;
}
} // namespace transition
