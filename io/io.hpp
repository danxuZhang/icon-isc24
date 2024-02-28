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
#include "../core/common/types.hpp"
#include <fstream>
#include <iostream>
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;

namespace io_muphys {

#ifdef __SINGLE_PRECISION
using NCreal_t = NcFloat;
#else
using NCreal_t = NcDouble;
#endif

void parse_args(std::string &file, size_t &itime, real_t &dt, real_t &qnc,
                int argc, char **argv);

void input_vector(NcFile &datafile, array_1d_t<real_t> &v,
                  const std::string input, size_t &ncells, size_t &nlev,
                  size_t itime);
void input_vector(NcFile &datafile, array_1d_t<real_t> &v,
                  const std::string input, size_t &ncells, size_t &nlev);

void output_vector(NcFile &datafile, array_1d_t<NcDim> &dims,
                   const std::string output, array_1d_t<real_t> &v,
                   size_t &ncells, size_t &nlev, int &deflate_level);
void output_vector(NcFile &datafile, array_1d_t<NcDim> &dims,
                   const std::string output, std::map<std::string, NcVarAtt>,
                   array_1d_t<real_t> &v, size_t &ncells, size_t &nlev,
                   int &deflate_level);

void read_fields(const std::string input_file, size_t &itime, size_t &ncells,
                 size_t &nlev, array_1d_t<real_t> &z, array_1d_t<real_t> &t,
                 array_1d_t<real_t> &p, array_1d_t<real_t> &rho,
                 array_1d_t<real_t> &qv, array_1d_t<real_t> &qc,
                 array_1d_t<real_t> &qi, array_1d_t<real_t> &qr,
                 array_1d_t<real_t> &qs, array_1d_t<real_t> &qg);
void write_fields(const string output_file, size_t &ncells, size_t &nlev,
                  array_1d_t<real_t> &t, array_1d_t<real_t> &qv,
                  array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
                  array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
                  array_1d_t<real_t> &qg);
void write_fields(const string output_file, const string input_file,
                  size_t &ncells, size_t &nlev, array_1d_t<real_t> &t,
                  array_1d_t<real_t> &qv, array_1d_t<real_t> &qc,
                  array_1d_t<real_t> &qi, array_1d_t<real_t> &qr,
                  array_1d_t<real_t> &qs, array_1d_t<real_t> &qg);
void log_time(int64_t value);
} // namespace io_muphys
