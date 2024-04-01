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

constexpr real_t a_ct = 1.72; //  (15/32)*(PI**0.5)*(EIR/RHOW)*V0R*AR**(1/8)
constexpr real_t b_ct = static_cast<real_t>(7.0) / static_cast<real_t>(8.0);
constexpr real_t c_agg_ct = 2.46;
constexpr real_t b_agg_ct = 0.94878;

/**
* @brief Calculates the aggregation of ice by graupel
*
* @param [in] rho Density (kg/m^3)
* @param [in] qr Rain specific mass (kg/kg)
* @param [in] qg Graupel specific mass (kg/kg)
* @param [in] qi Ice specific mass (kg/kg)
* @param [in] sticking_eff Sticking efficiency (dimensionless)
* @return The aggregation rate of ice by graupel (kg/kg/s)
*/
TARGET real_t ice_to_graupel(real_t rho, real_t qr, real_t qg, real_t qi,
                             real_t sticking_eff) {

  real_t result = 0.0;
  if (qi > graupel_ct::qmin) {
    if (qg > graupel_ct::qmin) {
      result = sticking_eff * qi * c_agg_ct * pow(rho * qg, b_agg_ct);
    }
    if (qr > graupel_ct::qmin) {
      result = result + a_ct * qi * pow(rho * qr, b_ct);
    }
  }

  return result;
}

} // namespace transition
