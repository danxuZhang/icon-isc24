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

constexpr real_t a_coop = 5.000;  ///< parameter in cooper fit
constexpr real_t b_coop = 0.304;  ///< parameter in cooper fit
constexpr real_t nimax = 250.e+3; ///< maximal number of ice crystals

/**
* @brief Calculates the ice crystal number concentration following the Cooper parameterization.
*
* This function calculates the ice crystal number concentration based on the Cooper
* parameterization, which depends on the ambient temperature and density. The number
* concentration is constrained to be less than or equal to nimax.
*
* @param [in] t Ambient temperature (K).
* @param [in] rho Ambient density (kg/m^3).
* @return The ice crystal number concentration (1/kg).
*/
TARGET real_t ice_number(real_t t, real_t rho) {
  return fmin(nimax, a_coop * exp(b_coop * (thermodyn::tmelt - t))) / rho;
}

} // namespace property
