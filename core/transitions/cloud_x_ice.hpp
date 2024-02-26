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

/**
 * @brief TODO
 *
 * @param [in] t Temperature
 * @param [in] qc Cloud specific mass
 * @param [in] qi Ice specific mass
 * @param [in] dt Time step
 * @return Homogeneous freezing rate
 */
TARGET real_t cloud_x_ice(real_t t, real_t qc, real_t qi, real_t dt) {
  real_t result = 0.0;

  if (qc > graupel_ct::qmin && t < graupel_ct::tfrz_hom)
    result = qc / dt;

  if (qi > graupel_ct::qmin && t > thermodyn::tmelt)
    result = -qi / dt;

  return result;
}
} // namespace transition
