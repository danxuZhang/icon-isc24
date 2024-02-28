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

constexpr real_t a_freez = 0.09;   // scale factor for freezing depression
constexpr real_t b_max_exp = 1.00; // maximum for exponential temperature factor
constexpr real_t eff_min = 0.075;  // minimum sticking efficiency
constexpr real_t eff_fac =
    3.5E-3; // Scaling factor [1/K] for cloud ice sticking efficiency
constexpr real_t tcrit =
    thermodyn::tmelt -
    85.; //   Temperature at which cloud ice autoconversion starts

/**
 * @brief sticking efficiency of ice
 *
 * @param [in] t Temperature
 * @return Ice sticking
 */
TARGET real_t ice_sticking(real_t t) {

  // per original code seems like aggregation is allowed even with no snow
  // present
  return fmax(
      fmax(fmin(exp(a_freez * (t - thermodyn::tmelt)), b_max_exp), eff_min),
      eff_fac * (t - tcrit));
}
} // namespace property
