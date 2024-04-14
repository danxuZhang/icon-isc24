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

struct t_qx_ptr {
  t_qx_ptr() : p(nullptr), x(nullptr), sz(0) {}
  t_qx_ptr(array_1d_t<real_t> &p_, array_1d_t<real_t> &x_) : p(p_.data()), x(x_.data()), sz(p_.size()) {}

  real_t* p;
  real_t* x;
  size_t sz;
}; // pointer vector

template<typename T>
inline auto get_min_element(T* arr, size_t sz) -> T{
  T min_ele = arr[0];
  for (size_t i = 0; i < sz; ++i) {
    if (arr[i] < min_ele) {
      min_ele = arr[i];
    }
  }
  return min_ele;
}

template<typename T>
inline auto get_max_5(T a, T b, T c, T d, T e) -> T {
  T max1 = std::max(a, b);
  T max2 = std::max(max1, c);
  T max3 = std::max(max2, d);
  return std::max(max3, e);
}

template<typename T>
inline auto get_max_3(T a, T b, T c) -> T {
  T max1 = std::max(a, b);
  return std::max(max1, c);
}

/**
 * @brief TODO
 *
 * @param precip time step for integration of microphysics  (s)
 * @param params fall speed parameters
 * @param zeta dt/(2dz)
 * @param vc state dependent fall speed correction
 * @param flx flux into cell from above
 * @param vt terminal velocity
 * @param q specific mass of hydrometeor
 * @param q_kp1 specific mass in next lower cell
 * @param rho density
 */
#pragma acc routine seq
void precip(const real_t (&params)[3], real_t (&precip)[3], real_t zeta,
            real_t vc, real_t flx, real_t vt, real_t q, real_t q_kp1,
            real_t rho) {
  real_t rho_x, flx_eff, flx_partial;
  rho_x = q * rho;
  flx_eff = (rho_x / zeta) + static_cast<real_t>(2.0) * flx;
  flx_partial = rho_x * vc * fall_speed(rho_x, params);
  flx_partial = std::fmin(flx_partial, flx_eff);
  precip[0] = (zeta * (flx_eff - flx_partial)) /
              ((static_cast<real_t>(1.0) + zeta * vt) * rho); // q update
  precip[1] =
      (precip[0] * rho * vt + flx_partial) * static_cast<real_t>(0.5); // flx
  rho_x = (precip[0] + q_kp1) * static_cast<real_t>(0.5) * rho;
  precip[2] = vc * fall_speed(rho_x, params); // vt
}

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

  auto is_sig_present = (bool*) malloc(nvec * ke * sizeof(bool)); // is snow, ice or graupel present?

  auto ind_k = (size_t*) malloc(nvec * ke * sizeof(size_t)); // k index of gathered point
  auto ind_i = (size_t*) malloc(nvec * ke * sizeof(size_t)); // iv index of gathered point

  auto kmin = (size_t**) malloc(nvec * sizeof(size_t*));
  for (size_t i = 0; i < nvec; ++i) {
    kmin[i] = (size_t*) malloc(np * sizeof(size_t));
  } // first level with condensate

  auto eflx = (real_t*) malloc(nvec * sizeof(real_t)); // internal energy flux from precipitation (W/m2 )
  real_t** vt = (real_t**) malloc(nvec * sizeof(real_t*));
  for (size_t i = 0; i < nvec; ++i) {
    vt[i] = (real_t*) malloc(np * sizeof(real_t));
    for (size_t j = 0; j < np; ++j) {
      vt[i][j] = 0.0;
    }
  }

  t_qx_ptr q[6]; // vector of pointers to point to four hydrometeor inouts
  array_1d_t<real_t> emptyArray;
  q[0] = {prr_gsp, qr};
  q[1] = {pri_gsp, qi};
  q[2] = {prs_gsp, qs};
  q[3] = {prg_gsp, qg};

  q[4] = {emptyArray, qc};
  q[5] = {emptyArray, qv};

  size_t jmx_ = 0;

  // The loop is intentionally i<nlev; since we are using an unsigned integer
  // data type, when i reaches 0, and you try to decrement further, (to -1), it
  // wraps to the maximum value representable by size_t.
  #pragma acc parallel vector_length(128) async
  #pragma acc loop seq 
  for  (size_t it = 0; it < ke; ++it) {
    const size_t i = ke - 1 - it;
    size_t jmx = 0;
    #pragma acc loop gang vector
    for (size_t j = ivstart; j < ivend; j++) {
      size_t oned_vec_index = i * ivend + j;
      if ((get_max_5(q[lqc].x[oned_vec_index], q[lqr].x[oned_vec_index],
                     q[lqs].x[oned_vec_index], q[lqi].x[oned_vec_index],
                     q[lqg].x[oned_vec_index]) > qmin) or
          ((t[oned_vec_index] < tfrz_het2) and
           (q[lqv].x[oned_vec_index] >
            qsat_ice_rho(t[oned_vec_index], rho[oned_vec_index])))) {
        #pragma acc atomic capture 
        {
        jmx = jmx_;
        jmx_ = jmx_ + 1;
        }
        ind_k[jmx] = i;
        ind_i[jmx] = j;
        is_sig_present[jmx] =
            get_max_3(q[lqs].x[oned_vec_index], q[lqi].x[oned_vec_index],
                      q[lqg].x[oned_vec_index]) > qmin;
      }

      for (size_t ix = 0; ix < np; ix++) {
        if (i == (ke - 1)) {
          kmin[j][qp_ind[ix]] = ke + 1;
          q[qp_ind[ix]].p[j] = 0.0;
          vt[j][ix] = 0.0;
        }

        if (q[qp_ind[ix]].x[oned_vec_index] > qmin) {
          kmin[j][qp_ind[ix]] = i;
        }
      }
    }
  }

  #pragma acc parallel vector_length(128) async loop gang vector 
  for (size_t j = 0; j < jmx_; j++) {
    real_t sx2x[nx][nx];
    for (size_t i = 0; i < nx; ++i) {
      for (size_t j_ = 0; j_ < nx; ++j_) {
        sx2x[i][j_] = 0.0;
      }
    }

    real_t sink[nx];
    real_t dqdt[nx];

    size_t k = ind_k[j];
    size_t iv = ind_i[j];
    size_t oned_vec_index = k * ivend + iv;

    real_t dvsw = q[lqv].x[oned_vec_index] -
           qsat_rho(t[oned_vec_index], rho[oned_vec_index]);
    real_t qvsi = qsat_ice_rho(t[oned_vec_index], rho[oned_vec_index]);
    real_t dvsi = q[lqv].x[oned_vec_index] - qvsi;
    real_t n_snow = snow_number(t[oned_vec_index], rho[oned_vec_index],
                         q[lqs].x[oned_vec_index]);
    real_t l_snow = snow_lambda(rho[oned_vec_index], q[lqs].x[oned_vec_index], n_snow);

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

    real_t n_ice, m_ice, x_ice;
    real_t eta, ice_dep;

    if (t[oned_vec_index] < tmelt) {
      n_ice = ice_number(t[oned_vec_index], rho[oned_vec_index]);
      m_ice = ice_mass(q[lqi].x[oned_vec_index], n_ice);
      x_ice = ice_sticking(t[oned_vec_index]);

      if (is_sig_present[j]) {
        eta = deposition_factor(
            t[oned_vec_index],
            qvsi); // neglect cloud depth cor. from gcsp_graupel
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
    } else {
      sx2x[lqc][lqr] = sx2x[lqc][lqr] + sx2x[lqc][lqs] + sx2x[lqc][lqg];
      sx2x[lqc][lqs] = 0.0;
      sx2x[lqc][lqg] = 0.0;
      ice_dep = 0.0;
      eta = 0.0;
    }

    if (is_sig_present[j]) {
      real_t dvsw0 = q[lqv].x[oned_vec_index] - qsat_rho(tmelt, rho[oned_vec_index]);
      sx2x[lqv][lqs] =
          vapor_x_snow(t[oned_vec_index], p[oned_vec_index],
                       rho[oned_vec_index], q[lqs].x[oned_vec_index], n_snow,
                       l_snow, eta, ice_dep, dvsw, dvsi, dvsw0, dt);
      sx2x[lqs][lqv] = -std::fmin(sx2x[lqv][lqs], 0.0);
      sx2x[lqv][lqs] = std::fmax(sx2x[lqv][lqs], 0.0);
      sx2x[lqv][lqg] = vapor_x_graupel(
          t[oned_vec_index], p[oned_vec_index], rho[oned_vec_index],
          q[lqg].x[oned_vec_index], dvsw, dvsi, dvsw0, dt);
      sx2x[lqg][lqv] = -std::fmin(sx2x[lqv][lqg], 0.0);
      sx2x[lqv][lqg] = std::fmax(sx2x[lqv][lqg], 0.0);
      sx2x[lqs][lqr] =
          snow_to_rain(t[oned_vec_index], p[oned_vec_index],
                       rho[oned_vec_index], dvsw0, q[lqs].x[oned_vec_index]);
      sx2x[lqg][lqr] =
          graupel_to_rain(t[oned_vec_index], p[oned_vec_index],
                          rho[oned_vec_index], dvsw0, q[lqg].x[oned_vec_index]);
    }

    for (size_t ix = 0; ix < nx; ix++) {
      sink[qx_ind[ix]] = 0.0;
      if ((is_sig_present[j]) or (qx_ind[ix] == lqc) or (qx_ind[ix] == lqv) or
          (qx_ind[ix] == lqr)) {
        for (size_t i = 0; i < nx; i++) {
          sink[qx_ind[ix]] = sink[qx_ind[ix]] + sx2x[qx_ind[ix]][i];
        }
        real_t stot = q[qx_ind[ix]].x[oned_vec_index] / dt;

        if ((sink[qx_ind[ix]] > stot) &&
            (q[qx_ind[ix]].x[oned_vec_index] > qmin)) {
          real_t nextSink = 0.0;
          for (size_t i = 0; i < nx; i++) {
            sx2x[qx_ind[ix]][i] = sx2x[qx_ind[ix]][i] * stot / sink[qx_ind[ix]];
            nextSink = nextSink + sx2x[qx_ind[ix]][i];
          }
          sink[qx_ind[ix]] = nextSink;
        }
      }
    }

    for (size_t ix = 0; ix < nx; ix++) {
      real_t sx2x_sum = 0;
      for (size_t i = 0; i < nx; i++) {
        sx2x_sum = sx2x_sum + sx2x[i][qx_ind[ix]];
      }
      dqdt[qx_ind[ix]] = sx2x_sum - sink[qx_ind[ix]];
      q[qx_ind[ix]].x[oned_vec_index] = std::fmax(
          0.0, q[qx_ind[ix]].x[oned_vec_index] + dqdt[qx_ind[ix]] * dt);
    }

    real_t qice = q[lqs].x[oned_vec_index] + q[lqi].x[oned_vec_index] +
           q[lqg].x[oned_vec_index];
    real_t qliq = q[lqc].x[oned_vec_index] + q[lqr].x[oned_vec_index];
    real_t qtot = q[lqv].x[oned_vec_index] + qice + qliq;
    real_t cv = cvd + (cvv - cvd) * qtot + (clw - cvv) * qliq +
         (ci - cvv) * qice; // qtot? or qv?
    t[oned_vec_index] =
        t[oned_vec_index] +
        dt *
            ((dqdt[lqc] + dqdt[lqr]) * (lvc - (clw - cvv) * t[oned_vec_index]) +
             (dqdt[lqi] + dqdt[lqs] + dqdt[lqg]) *
                 (lsc - (ci - cvv) * t[oned_vec_index])) /
            cv;
  }


  const size_t k_end = (lrain) ? ke : kstart - 1;
  #pragma acc parallel async vector_length(128)
  #pragma acc loop seq 
  for (size_t k = kstart; k < k_end; k++) {
    # pragma acc loop gang vector
    for (size_t iv = ivstart; iv < ivend; iv++) {
      size_t oned_vec_index = k * ivend + iv;
      if (k == kstart) {
        eflx[iv] = 0.0;
      }

      size_t kp1 = std::min(ke - 1, k + 1);
      if (k >= get_min_element(kmin[iv], np)) {
        real_t qliq = q[lqc].x[oned_vec_index] + q[lqr].x[oned_vec_index];
        real_t qice = q[lqs].x[oned_vec_index] + q[lqi].x[oned_vec_index] +
               q[lqg].x[oned_vec_index];

        real_t e_int =
            internal_energy(t[oned_vec_index], q[lqv].x[oned_vec_index], qliq,
                            qice, rho[oned_vec_index], dz[oned_vec_index]) +
            eflx[iv];
        real_t zeta = dt / (2.0 * dz[oned_vec_index]);
        real_t xrho = std::sqrt(rho_00 / rho[oned_vec_index]);

        for (size_t ix = 0; ix < np; ix++) {
          if (k >= kmin[iv][qp_ind[ix]]) {
            real_t update[3];
            real_t vc = vel_scale_factor(qp_ind[ix], xrho, rho[oned_vec_index],
                                  t[oned_vec_index],
                                  q[qp_ind[ix]].x[oned_vec_index]);
            precip(params[qp_ind[ix]], update, zeta, vc, q[qp_ind[ix]].p[iv],
                   vt[iv][ix], q[qp_ind[ix]].x[oned_vec_index],
                   q[qp_ind[ix]].x[kp1 * ivend + iv], rho[oned_vec_index]);
            q[qp_ind[ix]].x[oned_vec_index] = update[0];
            q[qp_ind[ix]].p[iv] = update[1];
            vt[iv][ix] = update[2];
          }
        }

        pflx[oned_vec_index] = q[lqs].p[iv] + q[lqi].p[iv] + q[lqg].p[iv];
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
        t[oned_vec_index] =
            T_from_internal_energy(e_int, q[lqv].x[oned_vec_index], qliq, qice,
                                   rho[oned_vec_index], dz[oned_vec_index]);
      }
    }
  }

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
