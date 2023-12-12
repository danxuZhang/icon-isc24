#pragma once

#include "../common/constants.hpp"
#include "../common/types.hpp"
#include <cmath>

namespace transition {

constexpr real_t a_rim = 4.43;
constexpr real_t b_rim = 0.94878;

/**
 * @brief TODO
 *
 * @param [in] t Temperature
 * @param [in] rho Ambient density
 * @param [in] qc Snow specific mass
 * @param [in] qg Graupel specific mass
 * @return Graupel riming rate
 */
TARGET real_t cloud_to_graupel(real_t t, real_t rho, real_t qc, real_t qg) {

  return (fmin(qc, qg) > graupel_ct::qmin && t > graupel_ct::tfrz_hom)
             ? a_rim * qc * pow(qg * rho, b_rim)
             : static_cast<real_t>(0.0);
}
} // namespace transition
