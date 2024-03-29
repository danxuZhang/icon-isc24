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


/**
 * @brief Test fixture class for Muphys tests.
 *
 * This class provides utility functions for validating floating-point values
 * and defines a constant value for zero.
 */
class MuphysTest : public testing::Test {

public:
      /**
     * @brief Validates the equality of actual and expected floating-point values.
     *
     * This function compares the actual and expected values based on the type of `real_t`.
     *
     * @param actual The actual value to be compared.
     * @param expected The expected value to compare against.
     */
  static void validate(real_t actual, real_t expected) {
    if constexpr (std::is_same_v<real_t, float>) {
      EXPECT_FLOAT_EQ(expected, actual);
    } else {
      EXPECT_DOUBLE_EQ(expected, actual);
    }
  }
  /**
     * @brief Validates the equality of actual and expected floating-point values.
     *
     * This function compares the actual value with the expected values based on the type of `real_t`.
     *
     * @param actual The actual value to be compared.
     * @param expected_float The expected value to compare against if `real_t` is `float`.
     * @param expected_double The expected value to compare against if `real_t` is `double`.
     */
  static void validate(real_t actual, real_t expected_float,
                       real_t expected_double) {
    if constexpr (std::is_same_v<real_t, float>) {
      EXPECT_FLOAT_EQ(expected_float, actual);
    } else {
      EXPECT_DOUBLE_EQ(expected_double, actual);
    }
  }
  /**
     * @brief A constant value representing zero.
     */
  static constexpr real_t ZERO = 0.0;
};
