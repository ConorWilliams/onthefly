#pragma once

#include "libatom/random/xoshiro.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::saddle {

  /**
   * @brief Provide a random pertubation to every atom within rcut of centre.
   *
   * Uses the minimum image distance to determine distance from centre. The pertubation is gaussian
   * in each coordinate axis and has standard deviation stddev.
   */
  void perturb(Vec3<floating> const& centre, SimCell& cell, floating rcut, floating stddev);

}  // namespace otf::saddle