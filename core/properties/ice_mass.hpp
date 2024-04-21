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

namespace property {

constexpr real_t mi_max = 1.0e-09; ///< maximum mass of cloud ice crystals (kg)

/**
* @brief Calculates the mass of an individual ice crystal.
*
* This function calculates the mass of an individual ice crystal based on the
* specific mass of ice and the ice crystal number concentration. The mass is
* constrained to be within a certain range defined by graupel_ct::m0_ice and mi_max.
*
* @param [in] qi Specific mass of ice (kg/kg).
* @param [in] ni Ice crystal number concentration (1/kg).
* @return The mass of an individual ice crystal (kg).
*/
TARGET real_t ice_mass(real_t qi, real_t ni) {
  return fmax(graupel_ct::m0_ice, fmin(qi / ni, mi_max));
}
} // namespace property
