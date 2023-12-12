#pragma once

#include "../common/constants.hpp"
#include "../common/types.hpp"
#include <cmath>

namespace property {

/**
 * @brief TODO
 * @param [in] density TODO
 * @param [in] params TODO
 * @return Fall speed
 */
template <typename array_t>
TARGET real_t fall_speed(real_t density, array_t params) {
  return params[0] * pow((density + params[2]), params[1]);
}
} // namespace property
