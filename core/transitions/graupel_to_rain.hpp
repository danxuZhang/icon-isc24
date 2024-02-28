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

#include "../common//types.hpp"
#include "../common/constants.hpp"
#include <cmath>

namespace transition {

constexpr real_t c1_melt = 12.31698;    // Constants in melting formula
constexpr real_t c2_melt = 7.39441e-05; // Constants in melting formula
constexpr real_t a_melt =
    graupel_ct::tx - static_cast<real_t>(389.5); // melting prefactor
constexpr real_t b_melt =
    static_cast<real_t>(3.0) / static_cast<real_t>(5.0); // melting exponent

/**
 * @brief Melting of graupel to form rain
 *
 * @param [in] t Ambient temperature
 * @param [in] p Ambient pressure
 * @param [in] rho Ambient density
 * @param [in] dvsw0 qv-qsat_water(T0)
 * @param [in] qg graupel specific mass
 * @return TODO
 */
TARGET real_t graupel_to_rain(real_t t, real_t p, real_t rho, real_t dvsw0,
                              real_t qg) {

  return (t > fmax(thermodyn::tmelt,
                   thermodyn::tmelt - graupel_ct::tx * dvsw0) &&
          qg > graupel_ct::qmin)
             ? (c1_melt / p + c2_melt) *
                   (t - thermodyn::tmelt + a_melt * dvsw0) *
                   pow(qg * rho, b_melt)
             : static_cast<real_t>(0.0);
}
} // namespace transition
