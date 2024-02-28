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

constexpr real_t mi_max = 1.0e-09; // maximum mass of cloud ice crystals

/**
 * @brief TODO
 * @param [in] qi Ice specific mass
 * @param [in] ni Ice crystal number
 * @return ice mass
 */
TARGET real_t ice_mass(real_t qi, real_t ni) {
  return fmax(graupel_ct::m0_ice, fmin(qi / ni, mi_max));
}
} // namespace property
