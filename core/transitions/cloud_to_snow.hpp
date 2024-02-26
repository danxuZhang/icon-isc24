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

constexpr real_t ecs =
    0.9; // Collection efficiency for snow collecting cloud water
constexpr real_t b_rim_ = -(graupel_ct::v1s + static_cast<real_t>(3.0));
constexpr real_t c_rim = static_cast<real_t>(2.61) * ecs *
                         graupel_ct::v0s; // (with pi*gam(v1s+3)/4 = 2.610)

/**
 * @brief TODO
 *
 * @param [in] t Temperature
 * @param [in] qc Cloud specific mass
 * @param [in] qs Snow specific mass
 * @param [in] ns Snow number
 * @param [in] lambda   Snow slope parameter (lambda)
 * @return Riming snow rate
 */
TARGET real_t cloud_to_snow(real_t t, real_t qc, real_t qs, real_t ns,
                            real_t lambda) {

  return (fmin(qc, qs) > graupel_ct::qmin && t > graupel_ct::tfrz_hom)
             ? (c_rim * ns) * qc * pow(lambda, b_rim_)
             : static_cast<real_t>(0.0);
}
} // namespace transition
