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
#include "core/common/types.hpp"
#include <gtest/gtest.h>

class MuphysTest : public testing::Test {

public:
  static void validate(real_t actual, real_t expected) {
    if constexpr (std::is_same_v<real_t, float>) {
      EXPECT_FLOAT_EQ(expected, actual);
    } else {
      EXPECT_DOUBLE_EQ(expected, actual);
    }
  }

  static void validate(real_t actual, real_t expected_float,
                       real_t expected_double) {
    if constexpr (std::is_same_v<real_t, float>) {
      EXPECT_FLOAT_EQ(expected_float, actual);
    } else {
      EXPECT_DOUBLE_EQ(expected_double, actual);
    }
  }

  static constexpr real_t ZERO = 0.0;
};
