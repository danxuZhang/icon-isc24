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

#include <cmath>

#include "../common/constants.hpp"
#include "../common/types.hpp"

using namespace thermodyn;

namespace thermo {

constexpr real_t c1es = 610.78;
constexpr real_t c2es = c1es * thermodyn::rd / thermodyn::rv;
constexpr real_t c3les = 17.269;
constexpr real_t c3ies = 21.875;
constexpr real_t c4les = 35.86;
constexpr real_t c4ies = 7.66;
constexpr real_t c5les = c3les * (thermodyn::tmelt - c4les);
constexpr real_t c5ies = c3ies * (thermodyn::tmelt - c4ies);

/**
* @brief Calculates the internal energy of a system.
*
* This function calculates the internal energy of a system based on the temperature,
* specific mass of water vapor, liquid phases, solid phases, density, and layer thickness.
* The internal energy is computed using the moist isometric specific heat and latent heats
* of vaporization and sublimation.
*
* @param [in] t Temperature (K).
* @param [in] qv Specific mass of water vapor (kg/kg).
* @param [in] qliq Specific mass of liquid phases (kg/kg).
* @param [in] qice Specific mass of solid phases (kg/kg).
* @param [in] rho Density (kg/m^3).
* @param [in] dz Layer thickness (m).
* @return The internal energy of the system (J/m^2).
*/
TARGET real_t internal_energy(real_t t, real_t qv, real_t qliq, real_t qice,
                              real_t rho, real_t dz) {

  // total water specific mass
  real_t qtot = qliq + qice + qv;

  // moist isometric specific heat
  real_t cv = cvd * (static_cast<real_t>(1.0) - qtot) + cvv * qv + clw * qliq +
              graupel_ct::ci * qice;

  return rho * dz * (cv * t - qliq * graupel_ct::lvc - qice * graupel_ct::lsc);
}

/**
* @brief Calculates the temperature from the internal energy of a system.
*
* This function calculates the temperature of a system based on the internal energy,
* specific mass of water vapor, liquid phases, solid phases, density, and layer thickness.
* The temperature is computed by solving the internal energy equation for temperature,
* considering the moist isometric specific heat and latent heats of vaporization and sublimation.
*
* @param [in] U Internal energy of the system (J/m^2).
* @param [in] qv Specific mass of water vapor (kg/kg).
* @param [in] qliq Specific mass of liquid phases (kg/kg).
* @param [in] qice Specific mass of solid phases (kg/kg).
* @param [in] rho Density (kg/m^3).
* @param [in] dz Layer thickness (m).
* @return The temperature of the system (K).
*/
TARGET real_t T_from_internal_energy(real_t U, real_t qv, real_t qliq,
                                     real_t qice, real_t rho, real_t dz) {

  // total water specific mass
  real_t qtot = qliq + qice + qv;

  // moist isometric specific heat
  real_t cv = (cvd * (static_cast<real_t>(1.0) - qtot) + cvv * qv + clw * qliq +
               graupel_ct::ci * qice) *
              rho * dz;

  return (U + rho * dz * (qliq * graupel_ct::lvc + qice * graupel_ct::lsc)) /
         cv;
}

/**
* @brief Calculates the specific humidity from vapor pressure and total pressure.
*
* This function calculates the specific humidity based on the vapor pressure and total pressure,
* using the gas constants for dry air and water vapor.
*
* @param [in] pvapor Vapor pressure (Pa).
* @param [in] ptotal Total pressure (Pa).
* @return The specific humidity (kg/kg).
*/
[[maybe_unused]] TARGET real_t specific_humidity(real_t pvapor, real_t ptotal) {

  real_t rdv = rd / rv;
  real_t o_m_rdv = static_cast<real_t>(1.) - rdv;

  return rdv * pvapor / (ptotal - o_m_rdv * pvapor);
}

/**
* @brief Calculates the saturation pressure over water.
*
* This function calculates the saturation pressure over water based on the temperature,
* using an empirical formula.
*
* @param [in] t Temperature (K).
* @return The saturation pressure over water (Pa).
*/
TARGET real_t sat_pres_water(real_t t) {
  return c1es * exp(c3les * (t - thermodyn::tmelt) / (t - c4les));
}

/**
 * @brief Calculates saturation pressure over ice
 *
 * @param [in] t Temperature (kelvin)
 * @return Saturation pressure
 */
TARGET real_t sat_pres_ice(real_t t) {
  return c1es * exp(c3ies * (t - thermodyn::tmelt) / (t - c4ies));
}

/**
 * @brief Calculates saturation vapor pressure (over liquid) at constant density
 *
 * @param [in] t Temperature (kelvin)
 * @param [in] rho Density
 * @return saturation pressure
 */
TARGET real_t qsat_rho(real_t t, real_t rho) {
  return sat_pres_water(t) / (rho * rv * t);
}

/**
 * @brief Saturation vapor pressure (over ice) at constant density
 *
 * @param [in] t Temperature (kelvin)
 * @param [in] rho Density
 * @return saturation pressure
 */
TARGET real_t qsat_ice_rho(real_t t, real_t rho) {
  return sat_pres_ice(t) / (rho * rv * t);
}

/**
 * @brief Computes the derivative d(qsat_rho)/dT
 *
 * @param [in] qs Saturation vapor pressure (over liquid)
 * @param [in] t Temperature (kelvin)
 * @return derivative d(qsat_rho)/dT
 */
[[maybe_unused]] TARGET real_t dqsatdT_rho(real_t qs, real_t t) {
  return qs * (c5les / pow(t - c4les, static_cast<real_t>(2.0)) -
               static_cast<real_t>(1.0) / t);
}

/**
* @brief Computes the derivative of saturation vapor pressure (over liquid) with respect to temperature.
*
* This function computes the derivative of saturation vapor pressure over liquid water with respect to temperature,
* based on the saturation vapor pressure and temperature.
*
* @param [in] qs Saturation vapor pressure over liquid water (Pa).
* @param [in] t Temperature (K).
* @return The derivative of saturation vapor pressure with respect to temperature (Pa/K).
*/
[[maybe_unused]] TARGET real_t dqsatdT(real_t qs, real_t t) {
  return c5les * (static_cast<real_t>(1.0) + vtmpc1 * qs) * qs /
         pow((t - c4les), static_cast<real_t>(2));
}

/**
* @brief Computes the derivative of saturation vapor pressure (over ice) with respect to temperature.
*
* This function computes the derivative of saturation vapor pressure over ice with respect to temperature,
* based on the saturation vapor pressure and temperature.
*
* @param [in] qs Saturation vapor pressure over ice (Pa).
* @param [in] t Temperature (K).
* @return The derivative of saturation vapor pressure over ice with respect to temperature (Pa/K).
*/
[[maybe_unused]] TARGET real_t dqsatdT_ice(real_t qs, real_t t) {
  return c5ies * (static_cast<real_t>(1.0) + vtmpc1 * qs) * qs /
         pow((t - c4ies), static_cast<real_t>(2));
}

/**
 * @brief Computes internal energy of vaporization
 *
 * @param [in] t Temperature (kelvin)
 * @return Energy of vaporization
 */
[[maybe_unused]] TARGET real_t vaporization_energy(real_t t) {
  return graupel_ct::lvc + (cvv - clw) * t;
}

/**
 * @brief Computes internal energy of sublimation
 *
 * @param [in] t Temperature (kelvin)
 * @return Energy of sublimation
 */
[[maybe_unused]] TARGET real_t sublimation_energy(real_t t) {
  return als + (cpv - graupel_ct::ci) * (t - thermodyn::tmelt) - rv * t;
}

} // namespace thermo
