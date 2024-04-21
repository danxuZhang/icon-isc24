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
#include "snow_number.hpp"
#include <cmath>

namespace property {

constexpr real_t b_i = static_cast<real_t>(2.0) / static_cast<real_t>(3.0);
constexpr real_t b_s = -static_cast<real_t>(1.0) / static_cast<real_t>(6.0);

/**
* @brief Calculates the velocity scale factor for different types of condensate.
*
* This function calculates the velocity scale factor for different types of condensate
* based on the condensate type, square root of the density ratio, density of condensate,
* temperature, and specific mass. The velocity scale factor is used to adjust the fall
* velocity of the condensate.
*
* @param [in] iqx Index of the condensate type (lqi for cloud ice, lqs for snow, etc.).
* @param [in] xrho Square root of the density ratio (sqrt(rho_00/rho)).
* @param [in] rho Density of the condensate (kg/m^3).
* @param [in] t Temperature (K).
* @param [in] qx Specific mass of the condensate (kg/kg).
* @return The velocity scale factor for the specified condensate type (dimensionless).
*/
 #pragma acc routine seq
TARGET real_t vel_scale_factor(int iqx, real_t xrho, real_t rho, real_t t,
                               real_t qx) {

  switch (iqx) {
  case idx::lqi:
    return pow(xrho, b_i);
  case idx::lqs:
    return xrho * pow(snow_number(t, rho, qx), b_s);
  default:
    return xrho;
  }
}
} // namespace property
