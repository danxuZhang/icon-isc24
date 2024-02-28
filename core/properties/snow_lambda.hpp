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
    static_cast<real_t>(2.0);     // (with ams*gam(bms+1.0_wp) where gam(3) = 2)
constexpr real_t lmd_0 = 1.0e+10; // no snow value of lambda
constexpr real_t bx =
    static_cast<real_t>(1.0) / (graupel_ct::bms + static_cast<real_t>(1.0));
constexpr real_t qsmin_ = 0.0e-6;

/**
 * @brief TODO
 * @param [in] rho Ambient density
 * @param [in] qs Snow specific mass
 * @param [in] ns Snow number
 * @return riming snow rate
 */
TARGET real_t snow_lambda(real_t rho, real_t qs, real_t ns) {

  return (qs > graupel_ct::qmin) ? pow((a2 * ns / ((qs + qsmin_) * rho)), bx)
                                 : lmd_0;
}
} // namespace property
