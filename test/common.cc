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
#include <gtest/gtest.h>
#include <iostream>
#include <type_traits>

#include "MuphysTest.cc"
#include "core/common/graupel.hpp"
#include "core/common/utils.hpp"

TEST(CommonTest, CommonTestSuite_CheckPrecision) {
#ifdef __SINGLE_PRECISION
  EXPECT_EQ(typeid(float), typeid(real_t));
#else
  EXPECT_EQ(typeid(double), typeid(real_t));
#endif
}

TEST_F(MuphysTest, PropertyTestSuite_DepAutoConversion_default) {
  real_t qi = 2.02422e-23;
  real_t m_ice = 1e-12;
  real_t ice_dep = -2.06276e-05;

  real_t result = property::deposition_auto_conversion(qi, m_ice, ice_dep);
  validate(result, ZERO);
}

TEST_F(MuphysTest, PropertyTestSuite_DepAutoConversion) {
  real_t qi = 2.02422e-2;
  real_t m_ice = 1e-12;
  real_t ice_dep = -2.06276e-05;
  real_t reference = 6.6430804299795412e-08;

  real_t result = property::deposition_auto_conversion(qi + graupel_ct::qmin,
                                                       m_ice, -ice_dep);
  validate(result, reference);
}

TEST_F(MuphysTest, PropertyTestSuite_IceDepositionNucleation_default) {
  real_t t = 272.731;
  real_t qc = 0.0;
  real_t qi = 2.02422e-23;
  real_t ni = 5.05089;
  real_t dvsi = -0.000618828;
  real_t dt = 30.0;

  real_t result = property::ice_deposition_nucleation(t, qc, qi, ni, dvsi, dt);
  validate(result, ZERO);
}

TEST_F(MuphysTest, PropertyTestSuite_IceDepositionNucleation) {
  real_t t = 160.9;
  real_t qc = 1.0e-2;
  real_t qi = 2.02422e-23;
  real_t ni = 5.05089;
  real_t dvsi = 0.0001;
  real_t dt = 30.0;
  real_t reference = 1.6836299999999999e-13;

  real_t result = property::ice_deposition_nucleation(t, qc, qi, ni, dvsi, dt);
  validate(result, reference);
}

TEST_F(MuphysTest, PropertyTestSuite_DepFactor) {
  real_t t = 272.731;
  real_t qvsi = 0.00416891;
  real_t reference = 1.3234329478493952e-05;

  real_t result = property::deposition_factor(t, qvsi);
  validate(result, reference);
}

TEST_F(MuphysTest, PropertyTestSuite_FallSpeed) {
  real_t reference = 0.67882452435647411;

  real_t result = property::fall_speed(ZERO, params[0]);
  validate(result, reference);
}

TEST_F(MuphysTest, PropertyTestSuite_IceSticking) {
  real_t reference = 1.0;

  real_t result = property::ice_sticking(thermodyn::tmelt);
  validate(result, reference);
}

TEST_F(MuphysTest, PropertyTestSuite_IceNumber) {
  real_t t = 272.731;
  real_t rho = 1.12442;
  real_t reference_float = 5.0508094;
  real_t reference_double = 5.0507995893464388;

  real_t result = property::ice_number(t, rho);
  validate(result, reference_float, reference_double);
}

TEST_F(MuphysTest, PropertyTestSuite_IceMass) {
  real_t qi = 2.02422e-23;
  real_t ni = 5.05089;
  real_t reference = 1e-12;

  real_t result = property::ice_mass(qi, ni);
  validate(result, reference);
}

TEST_F(MuphysTest, PropertyTestSuite_SnowLambda_default) {
  real_t rho = 1.12204;
  real_t qs = graupel_ct::qmin;
  real_t ns = 1.76669e+07;

  real_t result = property::snow_lambda(rho, qs, ns);
  validate(result, property::lmd_0);
}

TEST_F(MuphysTest, PropertyTestSuite_SnowLambda) {
  real_t rho = 1.12204;
  real_t qs = graupel_ct::qmin + static_cast<real_t>(10e-12);
  real_t ns = 1.76669e+07;
  real_t reference_float = 601168.3125;
  real_t reference_double = 601168.04842091922;

  real_t result = property::snow_lambda(rho, qs, ns);
  validate(result, reference_float, reference_double);
}

TEST_F(MuphysTest, PropertyTestSuite_SnowNumber_default) {
  real_t t = 276.302;
  real_t rho = 1.17797;
  real_t qs = 8.28451e-24;

  real_t result = property::snow_number(t, rho, qs);
  validate(result, property::n0s0);
}

TEST_F(MuphysTest, PropertyTestSuite_SnowNumber) {
  real_t t = 276.302;
  real_t rho = 1.17797;
  real_t qs = 8.28451e-4;
  real_t reference = 3813750;

  real_t result = property::snow_number(t, rho, qs);
  validate(result, reference);
}

TEST_F(MuphysTest, PropertyTestSuite_VelScaleFactor_default) {
  real_t xrho = 1.16163;
  real_t rho = 0.907829;
  real_t t = 257.501;
  real_t qx = 4.24365e-06;
  real_t reference = 1.16163;

  real_t result = property::vel_scale_factor(idx::lqg, xrho, rho, t, qx);
  validate(result, reference);
}

TEST_F(MuphysTest, PropertyTestSuite_VelScaleFactor_ice) {
  real_t xrho = 1.17873;
  real_t rho = 0.881673;
  real_t t = 257.458;
  real_t qx = 1.27526e-06;
  real_t reference = 1.1158596098981044;

  real_t result = property::vel_scale_factor(idx::lqi, xrho, rho, t, qx);
  validate(result, reference);
}

TEST_F(MuphysTest, PropertyTestSuite_VelScaleFactor_snow) {
  real_t xrho = 1.17787;
  real_t rho = 0.882961;
  real_t t = 257.101;
  real_t qx = 5.78761e-06;
  real_t reference = 0.06633230453931642;

  real_t result = property::vel_scale_factor(idx::lqs, xrho, rho, t, qx);
  validate(result, reference);
}

TEST_F(MuphysTest, ThermoTestSuite_InternalEnergy) {
  real_t t = 255.756;
  real_t qv = 0.00122576;
  real_t qliq = 1.63837e-20;
  real_t qice = 1.09462e-08;
  real_t rho = 0.83444;
  real_t dz = 249.569;
  real_t reference = 38265357.270336017;

  real_t result = thermo::internal_energy(t, qv, qliq, qice, rho, dz);
  validate(result, reference);
}

TEST_F(MuphysTest, ThermoTestSuite_T_InternalEnergy) {
  real_t u = 38265357.270336017;
  real_t qv = 0.00122576;
  real_t qliq = 1.63837e-20;
  real_t qice = 1.09462e-08;
  real_t rho = 0.83444;
  real_t dz = 249.569;
  real_t reference = 255.75599999999997;

  real_t result = thermo::T_from_internal_energy(u, qv, qliq, qice, rho, dz);
  validate(result, reference);
}

TEST_F(MuphysTest, ThermoTestSuite_SatPresWater) {
  real_t t = 281.787;
  real_t reference = 1120.1604149806028;

  real_t result = thermo::sat_pres_water(t);
  validate(result, reference);
}

TEST_F(MuphysTest, ThermoTestSuite_SatPresIce) {
  real_t t = 281.787;
  real_t reference_float = 1216.774;
  real_t reference_double = 1216.7746246067475;

  real_t result = thermo::sat_pres_ice(t);
  validate(result, reference_float, reference_double);
}

TEST_F(MuphysTest, ThermoTestSuite_QSatRho) {
  real_t t = 281.787;
  real_t rho = 1.24783;
  real_t reference_float = 0.0069027566;
  real_t reference_double = 0.0069027592942577506;

  real_t result = thermo::qsat_rho(t, rho);
  validate(result, reference_float, reference_double);
}

TEST_F(MuphysTest, ThermoTestSuite_QSatIceRho) {
  real_t t = 281.787;
  real_t rho = 1.24783;
  real_t reference_float = 0.0074981214;
  real_t reference_double = 0.0074981245870634101;

  real_t result = thermo::qsat_ice_rho(t, rho);
  validate(result, reference_float, reference_double);
}

TEST_F(MuphysTest, TransitionTestSuite_CloudToGraupel_default) {
  real_t t = 281.787;
  real_t rho = 1.24783;
  real_t qc = 0.0;
  real_t qg = 1.03636e-25;

  real_t result = transition::cloud_to_graupel(t, rho, qc, qg);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_CloudToGraupel) {
  real_t t = 256.983;
  real_t rho = 0.909677;
  real_t qc = 8.60101e-06;
  real_t qg = 4.11575e-06;
  real_t reference = 2.7054723496793982e-10;

  real_t result = transition::cloud_to_graupel(t, rho, qc, qg);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_CloudToRain_default) {
  real_t t = 281.787;
  real_t qc = 0.0;
  real_t qr = 5.2312e-07;
  real_t nc = 100;

  real_t result = transition::cloud_to_rain(t, qc, qr, nc);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_CloudToRain) {
  real_t t = 267.25;
  real_t qc = 5.52921e-05;
  real_t qr = 2.01511e-12;
  real_t nc = 100;
  real_t reference_float = 0.0045578829;
  real_t reference_double = 0.0045484481075162512;

  real_t result = transition::cloud_to_rain(t, qc, qr, nc);
  validate(result, reference_float, reference_double);
}

TEST_F(MuphysTest, TransitionTestSuite_CloudToSnow_default) {
  real_t t = 281.787;
  real_t qc = 0.0;
  real_t qs = 3.63983e-40;
  real_t ns = 800000;
  real_t lambda = 1e+10;

  real_t result = transition::cloud_to_snow(t, qc, qs, ns, lambda);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_CloudToSnow) {
  real_t t = 256.571;
  real_t qc = 3.31476e-05;
  real_t qs = 7.47365e-06;
  real_t ns = 3.37707e+07;
  real_t lambda = 8989.78;
  real_t reference = 9.5431874564438999e-10;

  real_t result = transition::cloud_to_snow(t, qc, qs, ns, lambda);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_CloudXIce_default) {
  real_t t = 256.835;
  real_t qc = 0.0;
  real_t qi = 4.50245e-07;
  real_t dt = 30;

  real_t result = transition::cloud_x_ice(t, qc, qi, dt);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_CloudXIce_i) {
  real_t t = thermodyn::tmelt + static_cast<real_t>(1.0);
  real_t qc = 0.0;
  real_t qi = 4.50245e-07;
  real_t dt = 30;
  real_t reference = -1.5008166666666666e-08;

  real_t result = transition::cloud_x_ice(t, qc, qi, dt);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_CloudXIce_c) {
  real_t t = 236;
  real_t qc = 0.000198441;
  real_t qi = 1.55543e-19;
  real_t dt = 30;
  real_t reference = 6.6146999999999995e-06;

  real_t result = transition::cloud_x_ice(t, qc, qi, dt);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_GraupelToRain_default) {
  real_t t = 280.156;
  real_t p = 98889.4;
  real_t rho = 1.22804;
  real_t dvsw0 = -0.00167867;
  real_t qg = 1.53968e-17;

  real_t result = transition::graupel_to_rain(t, p, rho, dvsw0, qg);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_GraupelToRain) {
  real_t t = 280.156;
  real_t p = 98889.4;
  real_t rho = 1.22804;
  real_t dvsw0 = -0.00167867;
  real_t qg = 1.53968e-15;
  real_t reference_float = 5.9748441e-13;
  real_t reference_double = 5.9748142538569357e-13;

  real_t result = transition::graupel_to_rain(t, p, rho, dvsw0, qg);
  validate(result, reference_float, reference_double);
}

TEST_F(MuphysTest, TransitionTestSuite_IceToGraupel_default) {
  real_t rho = 1.12442;
  real_t qr = 1.34006e-17;
  real_t qg = 1.22571e-13;
  real_t qi = 2.02422e-23;
  real_t sticking_eff = 0.962987;

  real_t result = transition::ice_to_graupel(rho, qr, qg, qi, sticking_eff);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_IceToGraupel_r) {
  real_t rho = 1.04848;
  real_t qr = 6.00408e-13;
  real_t qg = 1.19022e-18;
  real_t qi = 1.9584e-08;
  real_t sticking_eff = 0.518393;
  real_t reference = 7.1049436957697864e-19;

  real_t result = transition::ice_to_graupel(rho, qr, qg, qi, sticking_eff);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_IceToGraupel_g) {
  real_t rho = 1.04848;
  real_t qr = 6.00408e-16;
  real_t qg = 1.19022e-05;
  real_t qi = 1.9584e-08;
  real_t sticking_eff = 0.518393;
  real_t reference = 5.557203647060206e-13;

  real_t result = transition::ice_to_graupel(rho, qr, qg, qi, sticking_eff);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_IceToSnow_default) {
  real_t qi = 7.95122e-25;
  real_t ns = 2.23336e+07;
  real_t lambda = 61911.1;
  real_t sticking_eff = 0.241568;

  real_t result = transition::ice_to_snow(qi, ns, lambda, sticking_eff);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_IceToSnow) {
  real_t qi = 6.43223e-08;
  real_t ns = 1.93157e+07;
  real_t lambda = 10576.8;
  real_t sticking_eff = 0.511825;
  real_t reference = 3.3262745200740486e-11;

  real_t result = transition::ice_to_snow(qi, ns, lambda, sticking_eff);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_RainToGraupel_default) {
  real_t t = 272.731;
  real_t rho = 1.12442;
  real_t qc = 0.0;
  real_t qr = 1.34006e-17;
  real_t qi = 2.02422e-23;
  real_t qs = 1.02627e-19;
  real_t mi = 1e-12;
  real_t dvsw = -0.000635669;
  real_t dt = 30;

  real_t result =
      transition::rain_to_graupel(t, rho, qc, qr, qi, qs, mi, dvsw, dt);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_RainToGraupel_1) {
  real_t t = 258.542;
  real_t rho = 0.956089;
  real_t qc = 0.0;
  real_t qr = 3.01332e-11;
  real_t qi = 5.57166e-06;
  real_t qs = 3.55432e-05;
  real_t mi = 1e-09;
  real_t dvsw = 0.0;
  real_t dt = 30;
  real_t reference = 5.5463006215784835e-17;

  real_t result =
      transition::rain_to_graupel(t, rho, qc, qr, qi, qs, mi, dvsw, dt);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_RainToGraupel_2) {
  real_t t = 230.542;
  real_t rho = 0.956089;
  real_t qc = 8.6157e-05;
  real_t qr = 3.01332e-11;
  real_t qi = 5.57166e-06;
  real_t qs = 3.55432e-05;
  real_t mi = 1e-09;
  real_t dvsw = 0.0;
  real_t dt = 30;
  real_t reference = 1.0044953165173373e-12;

  real_t result =
      transition::rain_to_graupel(t, rho, qc, qr, qi, qs, mi, dvsw, dt);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_RainToGraupel_3) {
  real_t t = 258.542;
  real_t rho = 0.956089;
  real_t qc = 8.6157e-05;
  real_t qr = 3.01332e-11;
  real_t qi = 5.57166e-06;
  real_t qs = 3.55432e-05;
  real_t mi = 1e-09;
  real_t dvsw = 0.0;
  real_t dt = 30;
  real_t reference = 5.5316517337096312e-17;

  real_t result =
      transition::rain_to_graupel(t, rho, qc, qr, qi, qs, mi, dvsw, dt);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_RainToVapor_default) {
  real_t t = 258.542;
  real_t rho = 0.956089;
  real_t qc = 8.6157e-05;
  real_t qr = 3.01332e-11;
  real_t dvsw = 0.0;
  real_t dt = 30;

  real_t result = transition::rain_to_vapor(t, rho, qc, qr, dvsw, dt);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_RainToVapor) {
  real_t t = 258.542;
  real_t rho = 0.956089;
  real_t qc = 0.0;
  real_t qr = 3.01332e-11;
  real_t dvsw = -1e-10;
  real_t dt = 30;
  real_t reference_float = 2.8556715e-19;
  real_t reference_double = 2.8556697055499901e-19;

  real_t result = transition::rain_to_vapor(t, rho, qc, qr, dvsw, dt);
  validate(result, reference_float, reference_double);
}

TEST_F(MuphysTest, TransitionTestSuite_SnowToGraupel_default) {
  real_t t = 258.157;
  real_t rho = 0.93171;
  real_t qc = 0.0;
  real_t qs = 4.34854e-05;

  real_t result = transition::snow_to_graupel(t, rho, qc, qs);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_SnowToGraupel) {
  real_t t = 265.85;
  real_t rho = 1.04848;
  real_t qc = 7.02792e-05;
  real_t qs = 4.44664e-07;
  real_t reference = 6.2696154545048011e-10;

  real_t result = transition::snow_to_graupel(t, rho, qc, qs);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_SnowToRain_default) {
  real_t t = 265.83;
  real_t p = 80134.5;
  real_t rho = 1.04892;
  real_t dvsw0 = -0.00258631;
  real_t qs = 1.47687e-06;

  real_t result = transition::snow_to_rain(t, p, rho, dvsw0, qs);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_SnowToRain) {
  real_t t = 275.83;
  real_t p = 80134.5;
  real_t rho = 1.04892;
  real_t dvsw0 = 0.00258631;
  real_t qs = 1.47687e-06;
  real_t reference_float = 3.7268515e-07;
  real_t reference_double = 3.7268547760462804e-07;

  real_t result = transition::snow_to_rain(t, p, rho, dvsw0, qs);
  validate(result, reference_float, reference_double);
}

TEST_F(MuphysTest, TransitionTestSuite_VaporXGraupel_default) {
  real_t t = 278.026;
  real_t p = 95987.1;
  real_t rho = 1.20041;
  real_t qg = 2.05496e-16;
  real_t dvsw = -0.00234674;
  real_t dvsi = -0.00261576;
  real_t dvsw0 = -0.00076851;
  real_t dt = 30;

  real_t result =
      transition::vapor_x_graupel(t, p, rho, qg, dvsw, dvsi, dvsw0, dt);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_VaporXGraupel) {
  real_t t = 278.026;
  real_t p = 95987.1;
  real_t rho = 1.20041;
  real_t qg = 2.05496e-11;
  real_t dvsw = -0.00234674;
  real_t dvsi = -0.00261576;
  real_t dvsw0 = -0.00076851;
  real_t dt = 30;
  real_t reference = -6.8498666666666675e-13;

  real_t result =
      transition::vapor_x_graupel(t, p, rho, qg, dvsw, dvsi, dvsw0, dt);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_VaporXSnow_default) {
  real_t t = 278.748;
  real_t p = 95995.5;
  real_t rho = 1.19691;
  real_t qs = 1.25653e-20;
  real_t ns = 800000;
  real_t lambda = 1e+10;
  real_t eta = 0.0;
  real_t ice_dep = 0.0;
  real_t dvsw = -0.00196781;
  real_t dvsi = -0.00229367;
  real_t dvsw0 = -0.000110022;
  real_t dt = 30;

  real_t result = transition::vapor_x_snow(t, p, rho, qs, ns, lambda, eta,
                                           ice_dep, dvsw, dvsi, dvsw0, dt);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_VaporXSnow) {
  real_t t = 278.748;
  real_t p = 95995.5;
  real_t rho = 1.19691;
  real_t qs = 1.25653e-10;
  real_t ns = 800000;
  real_t lambda = 1e+10;
  real_t eta = 0.0;
  real_t ice_dep = 0.0;
  real_t dvsw = -0.00196781;
  real_t dvsi = -0.00229367;
  real_t dvsw0 = -0.000110022;
  real_t dt = 30;
  real_t reference = -8.6584296264775935e-13;

  real_t result = transition::vapor_x_snow(t, p, rho, qs, ns, lambda, eta,
                                           ice_dep, dvsw, dvsi, dvsw0, dt);
  validate(result, reference);
}

TEST_F(MuphysTest, TransitionTestSuite_VaporXIce_default) {
  real_t qi = 2.02422e-23;
  real_t mi = 1e-12;
  real_t eta = 1.32343e-05;
  real_t dvsi = -0.000618828;
  real_t dt = 30;

  real_t result = transition::vapor_x_ice(qi, mi, eta, dvsi, dt);
  validate(result, ZERO);
}

TEST_F(MuphysTest, TransitionTestSuite_VaporXIce) {
  real_t qi = 9.53048e-07;
  real_t mi = 1e-09;
  real_t eta = 1.90278e-05;
  real_t dvsi = 0.000120375;
  real_t dt = 30;
  real_t reference_float = 1.8469367e-09;
  real_t reference_double = 1.8469360555606007e-09;

  real_t result = transition::vapor_x_ice(qi, mi, eta, dvsi, dt);
  validate(result, reference_float, reference_double);
}

TEST(UtilsTestSuite, CalcDz) {
  array_1d_t<real_t> z = {1.0, -2.0, 1.0, 4.0, 5.0, 6.0};
  array_1d_t<real_t> dz;
  size_t ncells = 3;
  size_t nlev = 2;
  array_1d_t<real_t> expected = {-3, -7, -5, -3, -7, -5};

  utils_muphys::calc_dz(z, dz, ncells, nlev);

  EXPECT_EQ(dz.size(), ncells * nlev);
  for (size_t i = 0; i < dz.size(); i++)
    EXPECT_DOUBLE_EQ(dz[i], expected[i]);
}
