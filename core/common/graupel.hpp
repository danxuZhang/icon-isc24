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

#include "constants.hpp"
#include "types.hpp"

#include "../transitions/cloud_to_graupel.hpp"
#include "../transitions/cloud_to_rain.hpp"
#include "../transitions/cloud_to_snow.hpp"
#include "../transitions/cloud_x_ice.hpp"
#include "../transitions/graupel_to_rain.hpp"
#include "../transitions/ice_to_graupel.hpp"
#include "../transitions/ice_to_snow.hpp"
#include "../transitions/rain_to_graupel.hpp"
#include "../transitions/rain_to_vapor.hpp"
#include "../transitions/snow_to_graupel.hpp"
#include "../transitions/snow_to_rain.hpp"
#include "../transitions/vapor_x_graupel.hpp"
#include "../transitions/vapor_x_ice.hpp"
#include "../transitions/vapor_x_snow.hpp"

#include "../properties/deposition_auto_conversion.hpp"
#include "../properties/deposition_factor.hpp"
#include "../properties/fall_speed.hpp"
#include "../properties/ice_deposition.hpp"
#include "../properties/ice_mass.hpp"
#include "../properties/ice_number.hpp"
#include "../properties/ice_sticking.hpp"
#include "../properties/snow_lambda.hpp"
#include "../properties/snow_number.hpp"
#include "../properties/thermo.hpp"
#include "../properties/vel_scale_factor.hpp"

/**
 * @brief Simulates the microphysical processes involving graupel.
 *
 * This function simulates the microphysical processes involving graupel, including interactions
 * with other hydrometeor types such as cloud water, rain, snow, and ice. It calculates the
 * production and depletion rates of graupel and updates the specific masses of the hydrometeors
 * and the temperature.
 *
 * @param [in] nvec Number of horizontal points
 * @param [in] ke Number of grid points in vertical direction
 * @param [in] ivstart Start index for horizontal direction
 * @param [in] ivend End index for horizontal direction
 * @param [in] kstart Start index for vertical direction
 * @param [in] dt Time step for integration of microphysics (s)
 * @param [in] dz Layer thickness of full levels (m)
 * @param [inout] t Temperature in Kelvin
 * @param [in] rho Density of moist air (kg/m3)
 * @param [in] p Pressure (Pa)
 * @param [inout] qv Specific water vapor content (kg/kg)
 * @param [inout] qc Specific cloud water content (kg/kg)
 * @param [inout] qi Specific cloud ice content (kg/kg)
 * @param [inout] qr Specific rain content (kg/kg)
 * @param [inout] qs Specific snow content  kg/kg)
 * @param [inout] qg Specific graupel content (kg/kg)
 * @param [in] qnc Cloud number concentration
 * @param [out] prr_gsp Precipitation rate of rain, grid-scale (kg/(m2*s))
 * @param [out] pri_gsp Precipitation rate of ice, grid-scale (kg/(m2*s))
 * @param [out] prs_gsp Precipitation rate of snow, grid-scale (kg/(m2*s))
 * @param [out] prg_gsp Precipitation rate of graupel, grid-scale (kg/(m2*s))
 * @param [out] pflx Total precipitation flux
 *
 */
void graupel(size_t &nvec, size_t &ke, size_t &ivstart, size_t &ivend,
             size_t &kstart, real_t &dt, array_1d_t<real_t> &dz,
             array_1d_t<real_t> &t, array_1d_t<real_t> &rho,
             array_1d_t<real_t> &p, array_1d_t<real_t> &qv,
             array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
             array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
             array_1d_t<real_t> &qg, real_t &qnc, array_1d_t<real_t> &prr_gsp,
             array_1d_t<real_t> &pri_gsp, array_1d_t<real_t> &prs_gsp,
             array_1d_t<real_t> &prg_gsp, array_1d_t<real_t> &pflx);
