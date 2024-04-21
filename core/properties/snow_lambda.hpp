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

constexpr real_t a2 =
    graupel_ct::ams *
    static_cast<real_t>(2.0);     ///< (with ams*gam(bms+1.0_wp) where gam(3) = 2)
constexpr real_t lmd_0 = 1.0e+10; ///< no snow value of lambda
constexpr real_t bx = 
    static_cast<real_t>(1.0) / (graupel_ct::bms + static_cast<real_t>(1.0));  ///< Constant = 1.0 / (graupel_ct::bms + 1.0) 
constexpr real_t qsmin_ = 0.0e-6; ///< Minimum value of snow specific mass used in the calculation of snow lambda (kg/kg).

/**
* @brief Calculates the snow lambda parameter.
*
* This function calculates the snow lambda parameter based on the ambient density, snow specific mass,
* and snow number concentration. The lambda parameter is used in the size distribution of snow particles.
* If the snow specific mass is less than or equal to graupel_ct::qmin, a constant value of lmd_0 is returned.
*
* @param [in] rho Ambient density (kg/m^3).
* @param [in] qs Snow specific mass (kg/kg).
* @param [in] ns Snow number concentration (1/kg).
* @return The snow lambda parameter (1/m).
*/
TARGET real_t snow_lambda(real_t rho, real_t qs, real_t ns) {

  return (qs > graupel_ct::qmin) ? pow((a2 * ns / ((qs + qsmin_) * rho)), bx)
                                 : lmd_0;
}
} // namespace property
