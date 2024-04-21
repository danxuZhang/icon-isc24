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

constexpr real_t a_rim = 4.43;
constexpr real_t b_rim = 0.94878;

/**
 * @brief Calculates the graupel riming rate from cloud water.
 *
 * This function calculates the rate at which cloud water is converted to graupel through the
 * riming process. The riming rate is based on the temperature, ambient density, and the
 * specific masses of cloud water and graupel. If the minimum of cloud water and graupel
 * specific masses is below a threshold (graupel_ct::qmin) or the temperature is below the
 * homogeneous freezing threshold (graupel_ct::tfrz_hom), the riming rate is set to zero.
 *
 * @param [in] t Temperature (K).
 * @param [in] rho Ambient density (kg/m^3).
 * @param [in] qc Cloud water specific mass (kg/kg).
 * @param [in] qg Graupel specific mass (kg/kg).
 * @return The Graupel riming rate (kg/kg/s).
 */
TARGET real_t cloud_to_graupel(real_t t, real_t rho, real_t qc, real_t qg) {

  return (fmin(qc, qg) > graupel_ct::qmin && t > graupel_ct::tfrz_hom)
             ? a_rim * qc * pow(qg * rho, b_rim)
             : static_cast<real_t>(0.0);
}
} // namespace transition
