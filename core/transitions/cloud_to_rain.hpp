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

constexpr real_t qmin_ac = 1.00e-06; // threshold for auto conversion
constexpr real_t tau_max = 0.90e+00; // maximum allowed value of tau
constexpr real_t tau_min = 1.00e-30; // maximum allowed value of tau
constexpr real_t a_phi =
    6.00e+02; // constant in phi-function for autoconversion
constexpr real_t b_phi =
    0.68e+00; // exponent in phi-function for autoconversion
constexpr real_t c_phi = 5.00e-05;     // exponent in phi-function for accretion
constexpr real_t ac_kernel = 5.25e+00; // kernel coeff for SB2001 accretion
constexpr real_t x3 = 2.00e+00;        // gamma exponent for cloud distribution
constexpr real_t x2 = 2.60e-10;        // separating mass between cloud and rain
constexpr real_t x1 = 9.44e+09;        // kernel coeff for SB2001 autoconversion

/**
 * @brief TODO
 * @details
 * Kessler (1969) autoconversion rate
 *    scau = zccau * MAX( qc_ik - qc0, 0.0 )
 *    scac = zcac  * qc_ik * zeln7o8qrk
 * Seifert and Beheng (2001) autoconversion rate
 * with constant cloud droplet number concentration qnc
 *
 * @param [in] t Temperature
 * @param [in] qc Cloud water specific mass
 * @param [in] qr Rain water specific mass
 * @param [in] nc Cloud water number concentration
 * @return conversion rate
 */
TARGET real_t cloud_to_rain(real_t t, real_t qc, real_t qr, real_t nc) {
  const real_t au_kernel =
      x1 / (static_cast<real_t>(20.0) * x2) * (x3 + static_cast<real_t>(2.0)) *
      (x3 + static_cast<real_t>(4.0)) /
      pow((x3 + static_cast<real_t>(1.0)), static_cast<real_t>(2.0));
  real_t result = 0.0;

  if (qc > qmin_ac && t > graupel_ct::tfrz_hom) {
    real_t tau = fmax(tau_min, fmin(static_cast<real_t>(1.0) - qc / (qc + qr),
                                    tau_max)); // time-scale
    real_t phi = pow(tau, b_phi); // similarity function for autoconversion
    phi = a_phi * phi *
          pow((static_cast<real_t>(1.0) - phi), static_cast<real_t>(3.0));
    real_t xau = au_kernel * pow(qc * qc / nc, static_cast<real_t>(2.)) *
                 (static_cast<real_t>(1.0) +
                  phi / pow(static_cast<real_t>(1.0) - tau,
                            static_cast<real_t>(2.0))); // autoconversion rate
    real_t xac =
        ac_kernel * qc * qr *
        pow((tau / (tau + c_phi)), static_cast<real_t>(4.0)); // accretion rate
    result = xau + xac;
  }

  return result;
}

} // namespace transition
