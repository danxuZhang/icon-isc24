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

constexpr real_t a_rim_ct = 0.5; /// Constants in riming formula
constexpr real_t b_rim_ct = static_cast<real_t>(3.0) / static_cast<real_t>(4.0);

/**
* @brief Calculates the conversion rate from snow to graupel through riming
*
* @param [in] t Ambient temperature (K)
* @param [in] rho Ambient density (kg/m^3)
* @param [in] qc Cloud specific mass (kg/kg)
* @param [in] qs Snow specific mass (kg/kg)
* @returns The conversion rate from snow to graupel (kg/kg/s)
*/
TARGET real_t snow_to_graupel(real_t t, real_t rho, real_t qc, real_t qs) {

  return (fmin(qc, qs) > graupel_ct::qmin && t > graupel_ct::tfrz_hom)
             ? a_rim_ct * qc * pow(qs * rho, b_rim_ct)
             : static_cast<real_t>(0.0);
}

} // namespace transition
