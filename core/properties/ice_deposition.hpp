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
 * @brief Ice deposition nucleation
 * @param [in] t Temperature
 * @param [in] qc Specific humidity of ice
 * @param [in] qi Specific humidity of ice
 * @param [in] ni Ice crystal number
 * @param [in] dvsi Vapor excess with respect to ice sat
 * @param [in] dt Time step
 * @return Rate of vapor deposition for new ice
 */
TARGET real_t ice_deposition_nucleation(real_t t, real_t qc, real_t qi,
                                        real_t ni, real_t dvsi, real_t dt) {

  return (qi <= graupel_ct::qmin &&
          ((t < graupel_ct::tfrz_het2 && dvsi > 0.0) ||
           (t <= graupel_ct::tfrz_het1 && qc > graupel_ct::qmin)))
             ? fmin(graupel_ct::m0_ice * ni, fmax(0.0, dvsi)) / dt
             : static_cast<real_t>(0.0);
}

} // namespace property
