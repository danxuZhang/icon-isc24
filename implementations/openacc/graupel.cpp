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
// In numerical weather prediction (NWP) models like the one you're asking about, 
// the atmosphere is divided into a three-dimensional grid. Each grid point or "cell" 
// represents a discrete volume of air with its own set of atmospheric properties 
// (such as temperature, pressure, humidity) and, in this context, hydrometeor content. 
// Calculations are performed for each cell, considering interactions with neighboring cells, 
// to simulate the dynamics of the atmosphere over time. The term "cell" is commonly used 
// to denote these discrete elements of the computational grid.


#include "core/common/graupel.hpp"
#include "types.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <chrono>
#include <cmath>

#include <openacc.h>

using namespace property;
using namespace thermo;
using namespace transition;
using namespace idx;
using namespace graupel_ct;

/**
 * @struct t_qx_ptr
 * @brief Holds pointers to data arrays and their size for various hydrometeor types.
 * 
 * This structure is designed to manage arrays related to the properties of hydrometeors,
 * such as concentration, mass, volume, etc. It simplifies passing these arrays around
 * by encapsulating the pointers to the data and the size of the arrays in a single object.
 */
struct t_qx_ptr {
  real_t* p; ///< Pointer to an array holding properties of hydrometeors.
  real_t* x; ///< Pointer to an array holding another set of properties, possibly related to a different aspect of the hydrometeors.
  size_t sz; ///< Size of the arrays pointed to by `p` and `x`.

  /**
   * @brief Default constructor initializes pointers to `nullptr` and size to `0`.
   */
  t_qx_ptr() : p(nullptr), x(nullptr), sz(0) {}

  /**
   * @brief Constructor that initializes the structure with existing arrays.
   * @param p_ Reference to an array of real numbers for the first set of properties.
   * @param x_ Reference to an array of real numbers for the second set of properties.
   */
  t_qx_ptr(array_1d_t<real_t> &p_, array_1d_t<real_t> &x_) : p(p_.data()), x(x_.data()), sz(p_.size()) {}
};

template<typename T>
/**
 * @brief Finds the minimum element in an array.
 * 
 * This template function iterates over an array of any type that supports comparison
 * and identifies the minimum value contained within the array. It's a generic function
 * that can operate on arrays of integers, floats, or any other comparable types.
 * 
 * @param arr Pointer to the first element of the array to be searched.
 * @param sz The size of the array.
 * @return The minimum element found in the array.
 * @tparam T The type of the elements in the array. This type must support comparison operations.
 */
inline auto get_min_element(T* arr, size_t sz) -> T{
  T min_ele = arr[0];
  for (size_t i = 0; i < sz; ++i) {
    if (arr[i] < min_ele) {
      min_ele = arr[i];
    }
  }
  return min_ele;
}

#pragma acc routine seq
/**
 * @brief Simulates the update of precipitation in a grid cell over a time step.
 *
 * This function models the dynamics of precipitation within a grid cell, including
 * the accumulation or reduction of precipitation mass, the flux of precipitation
 * into the cell from the cell directly above, and the calculation of terminal
 * velocity of precipitation particles. It takes into account the physical properties
 * of the air and precipitation, such as air density and fall speed parameters, to
 * update the state of precipitation in the cell and its impact on adjacent cells.
 *
 * @param params An array of real_t containing parameters related to the fall speed of precipitation.
 * @param precip An array of real_t where the updated precipitation values will be stored.
 *               - precip[0] will contain the updated specific mass of hydrometeors (how much precipitation is present).
 *               - precip[1] will contain the flux of precipitation into the cell from above.
 *               - precip[2] will contain the updated terminal velocity of the precipitation particles.
 * @param zeta The ratio of the time step for integration (dt) to twice the vertical distance between grid cells (2dz).
 * @param vc The correction factor for the fall speed, which is state dependent (e.g., adjusted for local atmospheric conditions).
 * @param flx The flux of precipitation entering the cell from the cell directly above at the start of the time step.
 * @param vt The terminal velocity of precipitation particles (assumed constant for the time step).
 * @param q The specific mass of hydrometeors (precipitation) in the cell at the start of the time step.
 * @param q_kp1 The specific mass of hydrometeors in the cell directly below at the start of the time step.
 * @param rho The density of air in the cell.
 */
void precip(const real_t (&params)[3], real_t (&precip)[3], real_t zeta,
            real_t vc, real_t flx, real_t vt, real_t q, real_t q_kp1,
            real_t rho) {
  real_t rho_x, flx_eff, flx_partial;
  // Calculation of the effective density of precipitation in the cell
  rho_x = q * rho;
  // Effective flux of precipitation into the cell, considering both the density and the initial flux
  flx_eff = (rho_x / zeta) + static_cast<real_t>(2.0) * flx;
  // Partial flux of precipitation, adjusted for the fall speed correction factor
  flx_partial = rho_x * vc * fall_speed(rho_x, params);
  flx_partial = std::fmin(flx_partial, flx_eff);
  // Update the specific mass of hydrometeors in the cell
  precip[0] = (zeta * (flx_eff - flx_partial)) /
              ((static_cast<real_t>(1.0) + zeta * vt) * rho);
  // Update the flux of precipitation from the cell above
  precip[1] = (precip[0] * rho * vt + flx_partial) * static_cast<real_t>(0.5);
  // Recalculate the density of precipitation, considering the updates
  rho_x = (precip[0] + q_kp1) * static_cast<real_t>(0.5) * rho;
  // Update the terminal velocity of precipitation particles
  precip[2] = vc * fall_speed(rho_x, params);
}




/**
 * @brief Models the dynamics of graupel, snow, and ice precipitation within a model grid.
 *
 * This function simulates the microphysical processes affecting graupel, snow, and ice 
 * within a grid-based model of the atmosphere. It accounts for the formation, growth,
 * and fallout of these hydrometeors through the atmospheric column. The function updates the
 * state of each cell in terms of hydrometeor mass, flux, and related properties, considering
 * atmospheric conditions like temperature and density.
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
 * @details
 * The function operates by first identifying cells where significant hydrometeor processes are
 * expected to occur in loop 1, then applying microphysical transformations for significant cells
 * in loops 2 and finally updateshydrometeor states (such as mass and number concentration) in loop 3.
 * These updates are based on ambient atmospheric conditions
 * and the interactions between different hydrometeor types.
 */
void graupel(size_t &nvec, size_t &ke, size_t &ivstart, size_t &ivend,
             size_t &kstart, real_t &dt, array_1d_t<real_t> &dz,
             array_1d_t<real_t> &t, array_1d_t<real_t> &rho,
             array_1d_t<real_t> &p, array_1d_t<real_t> &qv,
             array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
             array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
             array_1d_t<real_t> &qg, real_t &qnc, array_1d_t<real_t> &prr_gsp,
             array_1d_t<real_t> &pri_gsp, array_1d_t<real_t> &prs_gsp,
             array_1d_t<real_t> &prg_gsp, array_1d_t<real_t> &pflx) {
  std::cout << "openacc graupel" << std::endl;


  // boolean array that will record whether snow, ice, or graupel is present in each model cell.
  auto is_sig_present = (bool*) malloc(nvec * ke * sizeof(bool)); 
  

  // arrays that store the vertical (k) and horizontal (iv) indices of model cells 
  auto ind_k = (size_t*) malloc(nvec * ke * sizeof(size_t)); // k index of gathered point
  auto ind_i = (size_t*) malloc(nvec * ke * sizeof(size_t)); // iv index of gathered point
  


  // 2D array that tracks the lowest vertical level (kmin) at which condensate 
  // (water in its liquid or solid form) is present for each vector. 
  auto kmin = (size_t**) malloc(nvec * sizeof(size_t*));
  for (size_t i = 0; i < nvec; ++i) {
    kmin[i] = (size_t*) malloc(np * sizeof(size_t));
  } 
  
  // Define variables for physical quantities and state parameters involved in precipitation processing.
  real_t cv, vc, eta, zeta, qvsi, qice, qliq, qtot, dvsw, dvsw0, dvsi, n_ice,
      m_ice, x_ice, n_snow, l_snow, ice_dep, e_int, stot, xrho;

  // scratch array with output from precipitation step
  real_t update[3]; 
  
  // internal energy flux (in Watts per square meter, W/m²) 
  // associated with precipitation within each vector.
  auto eflx = (real_t*) malloc(nvec * sizeof(real_t)); 


  // Allocate and initialize a 2D array to store the terminal velocity of precipitation particles in each model cell.
  real_t** vt = (real_t**) malloc(nvec * sizeof(real_t*));
  for (size_t i = 0; i < nvec; ++i) {
    vt[i] = (real_t*) malloc(np * sizeof(real_t));
    for (size_t j = 0; j < np; ++j) {
      vt[i][j] = 0.0;
    }
  }
  
  // Initialize t_qx_ptr objects to manage arrays of precipitation data.
  t_qx_ptr q[6];
  array_1d_t<real_t> emptyArray;
  q[0] = {prr_gsp, qr}; // Rain
  q[1] = {pri_gsp, qi}; // Ice
  q[2] = {prs_gsp, qs}; // Snow
  q[3] = {prg_gsp, qg}; // Graupel
  q[4] = {emptyArray, qc}; // Cloud water
  q[5] = {emptyArray, qv}; // Water vapor

  // Initializes variables to be used as indices in loops for processing the data. 
  size_t jmx = 0;
  size_t jmx_ = jmx;



  /**
  * @brief Loop 1
  * Traverses the simulation grid to identify cells where significant atmospheric processes are occurring.
  * by iterating over the vertical and horizontal dimensions of the simulation grid, then
  * identifying grid cells where significant atmospheric processes are expected to occur. The loop starts from
  * the top layer and moves downward, slicing horizontal layers vertically.
  */
  auto loop1_start = std::chrono::steady_clock::now();
  size_t oned_vec_index;
  #pragma acc parallel 
  #pragma acc loop seq 
  for  (size_t it = 0; it < ke; ++it) {
    const size_t i = ke - 1 - it; // start from top layer and go down

    #pragma acc loop private(oned_vec_index, jmx)
    for (size_t j = ivstart; j < ivend; j++) { // now slice horizontal layers vertically
      oned_vec_index = i * ivend + j;
      if ((std::max({q[lqc].x[oned_vec_index], q[lqr].x[oned_vec_index],
                     q[lqs].x[oned_vec_index], q[lqi].x[oned_vec_index],
                     q[lqg].x[oned_vec_index]}) > qmin) or
          ((t[oned_vec_index] < tfrz_het2) and
           (q[lqv].x[oned_vec_index] >
            qsat_ice_rho(t[oned_vec_index], rho[oned_vec_index])))) { // Determines if the current cell (i, j) is significant based on two criteria:
        #pragma acc atomic capture //record significant cell, label as sig present
        {
        jmx = jmx_;
        jmx_ = jmx_ + 1;
        }
        ind_k[jmx] = i;
        ind_i[jmx] = j;
        is_sig_present[jmx] =
            std::max({q[lqs].x[oned_vec_index], q[lqi].x[oned_vec_index],
                      q[lqg].x[oned_vec_index]}) > qmin;
      } 

      #pragma acc loop seq
      for (size_t ix = 0; ix < np; ix++) {
        if (i == (ke - 1)) { // initializes values to defaults,
          kmin[j][qp_ind[ix]] = ke + 1;
          q[qp_ind[ix]].p[j] = 0.0;
          vt[j][ix] = 0.0;
        }

        if (q[qp_ind[ix]].x[oned_vec_index] > qmin) {   
          kmin[j][qp_ind[ix]] = i; // topmost layer where condensate present
        }
      }
    }
  }
  auto loop1_end = std::chrono::steady_clock::now();
  auto loop1_duration = std::chrono::duration_cast<std::chrono::milliseconds>(loop1_end - loop1_start);
  std::cout << "time for loop one: " << loop1_duration.count() << "ms" << std::endl;


  /**
  * @brief Loop 2 
  * simulates the microphysical processes within the atmosphere at significant cells.By 
  * iterating over the significant cells identified in loop 1 and performs calculations
  * to model the microphysical processes occurring within those cells. It considers various
  * atmospheric variables and hydrometeor interactions to update the state of each significant cell.
  */
  auto loop2_start = std::chrono::steady_clock::now();
  size_t k, iv;
  real_t sx2x_sum;
  #pragma acc parallel loop private(k, iv, oned_vec_index, dvsw, qvsi, dvsi, n_snow, l_snow) \
                            private(n_ice, m_ice, x_ice, eta, ice_dep, dvsw0) \
                            private(stot, qice, qliq, qtot, cv, sx2x_sum) 
  for (size_t j = 0; j < jmx_; j++) {
    real_t sx2x[nx][nx]; //two-dimensional array used to store conversion rates 
                         //between different types at the point being processed.
    // sx2x[i][j] could represent the rate at which hydrometeor type i is converted into type j within a given time step in the simulation.
    for (size_t i = 0; i < nx; ++i) {
      for (size_t j_ = 0; j_ < nx; ++j_) {
        sx2x[i][j_] = 0.0;
      }
    }

    real_t sink[nx]; // track the cumulative effects of all sink processes of each hydrometeor type
    real_t dqdt[nx]; // rate of net decrease over time

    k = ind_k[j]; // vertical index of significant cell 
    iv = ind_i[j]; // horizontal index of significant cell 
    oned_vec_index = k * ivend + iv; // flattens the two-dimensional grid location (k, iv) 
    // into a one-dimensional array index for accessing atmospheric data arrays.

    dvsw = q[lqv].x[oned_vec_index] -
           qsat_rho(t[oned_vec_index], rho[oned_vec_index]);
    dvsi = q[lqv].x[oned_vec_index] - qvsi;
    // Differences between actual vapor content and saturation vapor content (drive phase changes)

    qvsi = qsat_ice_rho(t[oned_vec_index], rho[oned_vec_index]);
    n_snow = snow_number(t[oned_vec_index], rho[oned_vec_index],
                         q[lqs].x[oned_vec_index]);
    l_snow = snow_lambda(rho[oned_vec_index], q[lqs].x[oned_vec_index], n_snow);

    // rates of phase changes
    sx2x[lqc][lqr] = cloud_to_rain(t[oned_vec_index], q[lqc].x[oned_vec_index],
                                   q[lqr].x[oned_vec_index], qnc);
    sx2x[lqr][lqv] = rain_to_vapor(t[oned_vec_index], rho[oned_vec_index],
                                   q[lqc].x[oned_vec_index],
                                   q[lqr].x[oned_vec_index], dvsw, dt);
    sx2x[lqc][lqi] = cloud_x_ice(t[oned_vec_index], q[lqc].x[oned_vec_index],
                                 q[lqi].x[oned_vec_index], dt);
    sx2x[lqi][lqc] = -std::fmin(sx2x[lqc][lqi], 0.0);
    sx2x[lqc][lqi] = std::fmax(sx2x[lqc][lqi], 0.0);
    sx2x[lqc][lqs] = cloud_to_snow(t[oned_vec_index], q[lqc].x[oned_vec_index],
                                   q[lqs].x[oned_vec_index], n_snow, l_snow);
    sx2x[lqc][lqg] =
        cloud_to_graupel(t[oned_vec_index], rho[oned_vec_index],
                         q[lqc].x[oned_vec_index], q[lqg].x[oned_vec_index]);

    // whether processes like deposition (vapor to ice) or melting occur.
    if (t[oned_vec_index] < tmelt) {
      n_ice = ice_number(t[oned_vec_index], rho[oned_vec_index]);
      m_ice = ice_mass(q[lqi].x[oned_vec_index], n_ice);
      x_ice = ice_sticking(t[oned_vec_index]);

      // checks if significant hydrometeors (like snow, ice, or graupel) are present 
      if (is_sig_present[j]) {
        eta = deposition_factor( // how quickly vapor can turn into ice or snow.
            t[oned_vec_index],
            qvsi); // neglect cloud depth cor. from gcsp_graupel
        // Rate at which water vapor converts to ice
        sx2x[lqv][lqi] =
            vapor_x_ice(q[lqi].x[oned_vec_index], m_ice, eta, dvsi, dt);
        sx2x[lqi][lqv] = -std::fmin(sx2x[lqv][lqi], 0.0);
        sx2x[lqv][lqi] = std::fmax(sx2x[lqv][lqi], 0.0);
        ice_dep = std::fmin(sx2x[lqv][lqi], dvsi / dt);

        sx2x[lqi][lqs] = deposition_auto_conversion(q[lqi].x[oned_vec_index],
                                                    m_ice, ice_dep);
        sx2x[lqi][lqs] = sx2x[lqi][lqs] + ice_to_snow(q[lqi].x[oned_vec_index],
                                                      n_snow, l_snow, x_ice);
        sx2x[lqi][lqg] = ice_to_graupel(
            rho[oned_vec_index], q[lqr].x[oned_vec_index],
            q[lqg].x[oned_vec_index], q[lqi].x[oned_vec_index], x_ice);
        sx2x[lqs][lqg] =
            snow_to_graupel(t[oned_vec_index], rho[oned_vec_index],
                            q[lqc].x[oned_vec_index], q[lqs].x[oned_vec_index]);
        sx2x[lqr][lqg] = rain_to_graupel(
            t[oned_vec_index], rho[oned_vec_index], q[lqc].x[oned_vec_index],
            q[lqr].x[oned_vec_index], q[lqi].x[oned_vec_index],
            q[lqs].x[oned_vec_index], m_ice, dvsw, dt);
      }
      sx2x[lqv][lqi] =
          sx2x[lqv][lqi] +
          ice_deposition_nucleation(t[oned_vec_index], q[lqc].x[oned_vec_index],
                                    q[lqi].x[oned_vec_index], n_ice, dvsi, dt);
    } else { // > tmelt
      sx2x[lqc][lqr] = sx2x[lqc][lqr] + sx2x[lqc][lqs] + sx2x[lqc][lqg];
      sx2x[lqc][lqs] = 0.0;
      sx2x[lqc][lqg] = 0.0;
      ice_dep = 0.0;
      eta = 0.0;
    }

    if (is_sig_present[j]) {
      dvsw0 = q[lqv].x[oned_vec_index] - qsat_rho(tmelt, rho[oned_vec_index]);
      // difference between the actual water vapor content (q[lqv].x[oned_vec_index]) 
      // and the saturation vapor pressure at the melting point temperature (tmelt), 
      // for the given air density (rho). 
      // Indicates how far the current state is from saturation.
      sx2x[lqv][lqs] =
          vapor_x_snow(t[oned_vec_index], p[oned_vec_index],
                       rho[oned_vec_index], q[lqs].x[oned_vec_index], n_snow,
                       l_snow, eta, ice_dep, dvsw, dvsi, dvsw0, dt);
      // ensure that the calculated transformation rates do not fall into physically impossible values
      sx2x[lqs][lqv] = -std::fmin(sx2x[lqv][lqs], 0.0);
      sx2x[lqv][lqs] = std::fmax(sx2x[lqv][lqs], 0.0);
      sx2x[lqv][lqg] = vapor_x_graupel(
          t[oned_vec_index], p[oned_vec_index], rho[oned_vec_index],
          q[lqg].x[oned_vec_index], dvsw, dvsi, dvsw0, dt);
      
      // ensure that the calculated transformation rates do not fall into physically impossible values
      sx2x[lqg][lqv] = -std::fmin(sx2x[lqv][lqg], 0.0);
      sx2x[lqv][lqg] = std::fmax(sx2x[lqv][lqg], 0.0);
      sx2x[lqs][lqr] =
          snow_to_rain(t[oned_vec_index], p[oned_vec_index],
                       rho[oned_vec_index], dvsw0, q[lqs].x[oned_vec_index]);
      sx2x[lqg][lqr] =
          graupel_to_rain(t[oned_vec_index], p[oned_vec_index],
                          rho[oned_vec_index], dvsw0, q[lqg].x[oned_vec_index]);
    }
    // After calculating the exchange rates between hydrometeor types, 
    // the loop updates the concentrations (dqdt)

    #pragma acc loop seq
    for (size_t ix = 0; ix < nx; ix++) {
      sink[qx_ind[ix]] = 0.0;
      if ((is_sig_present[j]) or (qx_ind[ix] == lqc) or (qx_ind[ix] == lqv) or
          (qx_ind[ix] == lqr)) {
        #pragma acc loop seq
        for (size_t i = 0; i < nx; i++) {
          sink[qx_ind[ix]] = sink[qx_ind[ix]] + sx2x[qx_ind[ix]][i];
        }
        stot = q[qx_ind[ix]].x[oned_vec_index] / dt;
        // check if total sink exceeds the available amount (stot) of a hydrometeor.
        if ((sink[qx_ind[ix]] > stot) &&
            (q[qx_ind[ix]].x[oned_vec_index] > qmin)) {
          real_t nextSink = 0.0;
          #pragma acc loop seq
          for (size_t i = 0; i < nx; i++) {
            sx2x[qx_ind[ix]][i] = sx2x[qx_ind[ix]][i] * stot / sink[qx_ind[ix]];
            nextSink = nextSink + sx2x[qx_ind[ix]][i];
          }
          sink[qx_ind[ix]] = nextSink;
        }
      }
    }

    // ensures conservation by making sure that increases in one hydrometeor type 
    // are balanced by decreases in others (e.g., cloud water, if rain is formed by condensation).
    #pragma acc loop seq
    for (size_t ix = 0; ix < nx; ix++) {
      sx2x_sum = 0;
      #pragma acc loop seq
      for (size_t i = 0; i < nx; i++) {
        sx2x_sum = sx2x_sum + sx2x[i][qx_ind[ix]];
      }
      dqdt[qx_ind[ix]] = sx2x_sum - sink[qx_ind[ix]];
      q[qx_ind[ix]].x[oned_vec_index] = std::fmax(
          0.0, q[qx_ind[ix]].x[oned_vec_index] + dqdt[qx_ind[ix]] * dt);
    }

    qice = q[lqs].x[oned_vec_index] + q[lqi].x[oned_vec_index] +
           q[lqg].x[oned_vec_index];
    qliq = q[lqc].x[oned_vec_index] + q[lqr].x[oned_vec_index];
    qtot = q[lqv].x[oned_vec_index] + qice + qliq;
    cv = cvd + (cvv - cvd) * qtot + (clw - cvv) * qliq +
         (ci - cvv) * qice; // qtot? or qv?
    // adjustment of temperature based on the latent heat release or absorption during phase changes
    t[oned_vec_index] =
        t[oned_vec_index] +
        dt *
            ((dqdt[lqc] + dqdt[lqr]) * (lvc - (clw - cvv) * t[oned_vec_index]) +
             (dqdt[lqi] + dqdt[lqs] + dqdt[lqg]) *
                 (lsc - (ci - cvv) * t[oned_vec_index])) /
            cv;
  }
  // basically update to the overall atmospheric state at each significant point, 
  auto loop2_end = std::chrono::steady_clock::now();
  auto loop2_duration = std::chrono::duration_cast<std::chrono::milliseconds>(loop2_end - loop2_start);
  std::cout << "time for loop two: " << loop2_duration.count() << "ms" << std::endl;

  size_t kp1;
  size_t k_end = (lrain) ? ke : kstart - 1;
  auto loop3_start = std::chrono::steady_clock::now();


  /**
  * @brief Loop 3
  * Simulates vertical transport and captures feedback between precipitation and atmospheric conditions.
  * by traversing each vertical and horizontal cell of the simulation grid, simulating the vertical
  * transport of precipitation and capturing the feedback between precipitation-phase changes and atmospheric
  * conditions. It adjusts hydrometeor concentrations and velocities, and updates temperature based on latent
  * heat exchange. 
  */
  #pragma acc parallel 
  #pragma acc loop seq 
  for (size_t k = kstart; k < k_end; k++) {
    # pragma acc loop private(oned_vec_index, kp1, qliq, qice, e_int, zeta, xrho, vc, update)
    for (size_t iv = ivstart; iv < ivend; iv++) {
      oned_vec_index = k * ivend + iv;
      if (k == kstart) {
        eflx[iv] = 0.0;
      }

      kp1 = std::min(ke - 1, k + 1);
      
      // feedback between precipitation processes and atmospheric conditions
      if (k >= get_min_element(kmin[iv], np)) {
        // calculates the liquid and ice water content at each point
        qliq = q[lqc].x[oned_vec_index] + q[lqr].x[oned_vec_index];
        qice = q[lqs].x[oned_vec_index] + q[lqi].x[oned_vec_index] +
               q[lqg].x[oned_vec_index];
        // Internal energy (e_int) based on temperature, hydrometeor content, and air density.
        e_int =
            internal_energy(t[oned_vec_index], q[lqv].x[oned_vec_index], qliq,
                            qice, rho[oned_vec_index], dz[oned_vec_index]) +
            eflx[iv];

        zeta = dt / (2.0 * dz[oned_vec_index]);
        xrho = std::sqrt(rho_00 / rho[oned_vec_index]);

        
        // simulate how precipitation and its associated energy 
        // move vertically through the atmosphere.
        for (size_t ix = 0; ix < np; ix++) {
          if (k >= kmin[iv][qp_ind[ix]]) {
            vc = vel_scale_factor(qp_ind[ix], xrho, rho[oned_vec_index],
                                  t[oned_vec_index],
                                  q[qp_ind[ix]].x[oned_vec_index]);

            // simulates precipitation processes, adjusting the vertical distribution of different hydrometeors
            precip(params[qp_ind[ix]], update, zeta, vc, q[qp_ind[ix]].p[iv],
                   vt[iv][ix], q[qp_ind[ix]].x[oned_vec_index],
                   q[qp_ind[ix]].x[kp1 * ivend + iv], rho[oned_vec_index]);
            q[qp_ind[ix]].x[oned_vec_index] = update[0];
            q[qp_ind[ix]].p[iv] = update[1];
            vt[iv][ix] = update[2];
          }
        }

        pflx[oned_vec_index] = q[lqs].p[iv] + q[lqi].p[iv] + q[lqg].p[iv];

        // updates the flux of energy (eflx) associated with precipitation, 
        eflx[iv] =
            dt * (q[lqr].p[iv] * (clw * t[oned_vec_index] -
                                  cvd * t[kp1 * ivend + iv] - lvc) +
                  pflx[oned_vec_index] * (ci * t[oned_vec_index] -
                                          cvd * t[kp1 * ivend + iv] - lsc));
        pflx[oned_vec_index] = pflx[oned_vec_index] + q[lqr].p[iv];
        qliq = q[lqc].x[oned_vec_index] + q[lqr].x[oned_vec_index];
        qice = q[lqs].x[oned_vec_index] + q[lqi].x[oned_vec_index] +
               q[lqg].x[oned_vec_index];
        e_int = e_int - eflx[iv];

        // The temperature at each grid point is adjusted based on the updated internal energy
        t[oned_vec_index] =
            T_from_internal_energy(e_int, q[lqv].x[oned_vec_index], qliq, qice,
                                   rho[oned_vec_index], dz[oned_vec_index]);
      }
    }
  }

  auto loop3_end = std::chrono::steady_clock::now();
  auto loop3_duration = std::chrono::duration_cast<std::chrono::milliseconds>(loop3_end - loop3_start);
  std::cout << "time for loop three: " << loop3_duration.count() << "ms" << std::endl;

  for (size_t i = 0; i < nvec; ++i) {
    free(vt[i]);
    free(kmin[i]);
  }
  free(is_sig_present);
  free(vt);
  free(kmin);
  free(ind_k);
  free(ind_i);
  free(eflx);
}
