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

constexpr real_t m0_s = 3.0e-9; // initial mass of snow crystals
constexpr real_t b_dep = static_cast<real_t>(2.0) / static_cast<real_t>(3.0);
constexpr real_t xcrit = 1.0; // threshold parameter

/**
 * @brief TODO
 * @param [in] qi Ice specific mass
 * @param [in] m_ice Ice crystal mass
 * @param [in] ice_dep Rate of ice deposition (some to snow)
 * @return TODO
 */
TARGET real_t deposition_auto_conversion(real_t qi, real_t m_ice,
                                         real_t ice_dep) {

  real_t result = 0.0;
  if (qi > graupel_ct::qmin) {
    real_t tau_inv = b_dep / (pow((m0_s / m_ice), b_dep) - xcrit);
    result = fmax(static_cast<real_t>(0.0), ice_dep) * tau_inv;
  }

  return result;
}
} // namespace property
