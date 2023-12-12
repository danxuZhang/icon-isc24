#pragma once

#include "types.hpp"
#include <cmath>

namespace thermodyn {

// Thermodynamic constants for the dry and moist atmosphere

// Dry air
constexpr real_t rd = 287.04; // [J/K/kg] gas constexprant
constexpr real_t cpd =
    1004.64; // [J/K/kg] specific heat at constexprant pressure
constexpr real_t cvd =
    cpd - rd; // [J/K/kg] specific heat at constexprant volume
constexpr real_t con_m = 1.50E-5;  // [m^2/s]  kinematic viscosity of dry air
constexpr real_t con_h = 2.20E-5;  // [m^2/s]  scalar conductivity of dry air
constexpr real_t con0_h = 2.40e-2; // [J/m/s/K]thermal conductivity of dry air
constexpr real_t eta0d = 1.717e-5; // [N*s/m2] dyn viscosity of dry air at tmelt

// H2O
// gas
constexpr real_t rv = 461.51; // [J/K/kg] gas constexprant for water vapor
constexpr real_t cpv =
    1869.46; // [J/K/kg] specific heat at constexprant pressure
constexpr real_t cvv =
    cpv - rv; // [J/K/kg] specific heat at constexprant volume
constexpr real_t dv0 =
    2.22e-5; // [m^2/s]  diff coeff of H2O vapor in dry air at tmelt
// liquid / water
constexpr real_t rhoh2o = 1000.; // [kg/m3]  density of liquid water
// solid / ice
constexpr real_t rhoice = 916.7; // [kg/m3]  density of pure ice

constexpr real_t cv_i = 2000.0;

// phase changes
constexpr real_t alv = 2.5008e6;  // [J/kg]   latent heat for vaporisation
constexpr real_t als = 2.8345e6;  // [J/kg]   latent heat for sublimation
constexpr real_t alf = als - alv; // [J/kg]   latent heat for fusion
constexpr real_t tmelt = 273.15;  // [K]      melting temperature of ice/snow
constexpr real_t t3 = 273.16;     // [K]      Triple point of water at 611hPa

// Auxiliary constexprants
constexpr real_t rdv = rd / rv;                                // [ ]
constexpr real_t vtmpc1 = rv / rd - static_cast<real_t>(1.);   // [ ]
constexpr real_t vtmpc2 = cpv / cpd - static_cast<real_t>(1.); // [ ]
constexpr real_t rcpv = cpd / cpv - static_cast<real_t>(1.);   // [ ]
constexpr real_t alvdcp = alv / cpd;                           // [K]
constexpr real_t alsdcp = als / cpd;                           // [K]
constexpr real_t rcpd = static_cast<real_t>(1.) / cpd;         // [K*kg/J]
constexpr real_t rcvd = static_cast<real_t>(1.) / cvd;         // [K*kg/J]
constexpr real_t rcpl = 3.1733; // cp_d / cp_l - 1

constexpr real_t clw = (rcpl + static_cast<real_t>(1.0)) *
                       cpd; // specific heat capacity of liquid water
constexpr real_t cv_v = (rcpv + static_cast<real_t>(1.0)) * cpd - rv;
} // namespace thermodyn

namespace graupel_ct {

constexpr real_t rho_00 = 1.225; // reference air density
constexpr real_t q1 = 8.e-6;
constexpr real_t qmin = 1.0e-15; // threshold for computation
constexpr real_t ams =
    0.069; // Formfactor in the mass-size relation of snow particles
constexpr real_t bms =
    2.0; // Exponent in the mass-size relation of snow particles
constexpr real_t v0s = 25.0; // prefactor in snow fall speed
constexpr real_t v1s = 0.5;  // Exponent in the terminal velocity for snow
constexpr real_t m0_ice =
    1.0e-12;                 // initial crystal mass for cloud ice nucleation
constexpr real_t ci = 2108.; // specific heat of ice
constexpr real_t tx = 3339.5;
constexpr real_t tfrz_het1 =
    thermodyn::tmelt -
    static_cast<real_t>(
        6.0); // temperature for het. freezing of cloud water with supersat
constexpr real_t tfrz_het2 =
    thermodyn::tmelt -
    static_cast<real_t>(25.0); // temperature for het. freezing of cloud water
constexpr real_t tfrz_hom =
    thermodyn::tmelt -
    static_cast<real_t>(37.0); // temperature for hom. freezing of cloud water
constexpr real_t lvc =
    thermodyn::alv -
    (thermodyn::cpv - thermodyn::clw) *
        thermodyn::tmelt; // invariant part of vaporization enthalpy
constexpr real_t lsc =
    thermodyn::als -
    (thermodyn::cpv - ci) *
        thermodyn::tmelt; // invariant part of vaporization enthalpy
} // namespace graupel_ct

namespace idx {
constexpr size_t nx = 6;  // number of water species
constexpr size_t np = 4;  // number of precipitating water species
constexpr size_t lqr = 0; // index for rain
constexpr size_t lqi = 1; // index for ice
constexpr size_t lqs = 2; // index for snow
constexpr size_t lqg = 3; // index for graupel
constexpr size_t lqc = 4; // index for cloud
constexpr size_t lqv = 5; // index for vapor

constexpr size_t qx_ind[] = {lqv, lqc, lqr, lqs, lqi, lqg};
constexpr size_t qp_ind[] = {lqr, lqi, lqs, lqg};
} // namespace idx

constexpr bool lrain = true; // switch for disabling rain
constexpr bool lcold = true; // switch for disabling freezing processes

constexpr real_t params[4][3] = {
    {14.58, 0.111, 1.0e-12},
    {1.25, 0.160, 1.0e-12},
    {57.80, static_cast<real_t>(0.5) / static_cast<real_t>(3.0), 1.0e-12},
    {12.24, 0.217, 1.0e-08}};
