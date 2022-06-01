#pragma once

#include "libatom/random/xoshiro.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::saddle {

  void perturb(Vec3<floating> const& centre, SimCell& cell, floating rcut, floating stddev);

}  // namespace otf::saddle