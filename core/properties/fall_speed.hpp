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

/**
* @brief Calculates the fall speed of hydrometeors.
*
* This function calculates the fall speed of hydrometeors based on their density
* and a set of parameters. The fall speed is determined using a power law relationship.
*
* @param [in] density The density of the hydrometeors (kg/m^3).
* @param [in] params An array of parameters used in the fall speed calculation.
*                    - params[0]: Coefficient of the power law relationship.
*                    - params[1]: Exponent of the power law relationship.
*                    - params[2]: Offset added to the density.
* @return The fall speed of the hydrometeors (m/s).
*/
template <typename array_t>
TARGET real_t fall_speed(real_t density, array_t params) {
  return params[0] * pow((density + params[2]), params[1]);
}
} // namespace property
