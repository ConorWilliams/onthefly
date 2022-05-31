#pragma once



#include "libatom/sim_cell.hpp"
#include "libatom/random/xoshiro.hpp"
#include "libatom/utils.hpp"

namespace otf::saddle {

   void perturb(std::size_t i, SimCell &, floating dr);

}  // namespace otf::saddle