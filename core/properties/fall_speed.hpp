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

/**
 * @brief TODO
 * @param [in] density TODO
 * @param [in] params TODO
 * @return Fall speed
 */
template <typename array_t>
TARGET real_t fall_speed(real_t density, array_t params) {
  return params[0] * pow((density + params[2]), params[1]);
}
} // namespace property
