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
#include <exception>
#include <string>
#include <vector>
using namespace std;

//  using Integer_4 = int;
//  using Integer_8 = long long int;

#ifdef __SINGLE_PRECISION
using real_t = float;
#else
using real_t = double;
#endif

template <class T> T type_converter(char *str_ptr, char **end) {
  if constexpr (is_same_v<float, T>) {
    return strtof(str_ptr, end);
  } else if constexpr (is_same_v<double, T>) {
    return strtod(str_ptr, end);
  } else if constexpr (is_same_v<long double, T>) {
    return strtold(str_ptr, end);
  }
  // induce compile time error if none of the cases are fulfilled
  else
    static_assert(not is_same_v<T, T>, "Invalid template type");
}

template <typename T> using array_1d_t = std::vector<T, allocator<T>>;

template <typename T>
using array_2d_t = std::vector<std::vector<T, allocator<T>>>;

#ifdef MU_DEVICE
#define TARGET __host__ __device__
#else
#define TARGET inline
#endif
