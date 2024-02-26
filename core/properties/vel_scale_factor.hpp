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
 * @brief TODO
 * @param [in] xrho TODO sqrt(rho_00/rho)
 * @param [in] rho Density of condensate
 * @param [in] t Temperature
 * @param [in] qx Specific mass
 * @return Scale factor
 */
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
