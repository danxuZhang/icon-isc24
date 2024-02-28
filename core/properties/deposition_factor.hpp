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

constexpr real_t kappa = 2.40e-2; // thermal conductivity of dry air
constexpr real_t b = 1.94;
constexpr real_t a = thermodyn::als * thermodyn::als / (kappa * thermodyn::rv);

/**
 * @brief TODO
 * @param [in] t Temperature
 * @param [in] qvsi Saturation (ice) specific vapor mass
 * @return Deposition factor
 */
TARGET real_t deposition_factor(real_t t, real_t qvsi) {
  real_t cx = static_cast<real_t>(2.22e-5) * pow(thermodyn::tmelt, (-b)) *
              static_cast<real_t>(101325.0);

  real_t x = cx / thermodyn::rd * pow(t, b - static_cast<real_t>(1.0));
  return x / (static_cast<real_t>(1.0) + a * x * qvsi / (t * t));
}
} // namespace property
