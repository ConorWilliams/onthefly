#pragma once

#include <cstddef>

#include "libatom/asserts.hpp"
#include "libatom/neighbour/neighbour_list.hpp"
#include "libatom/system/atom_array.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf {

  class EAM {
  public:
    // Get this potentials cut-off radius
    double rcut() const;

    // Compute energy
    template <typename... Mems> double energy(SimCell const &, NeighbourList const &);

    // Compute gradient, force on frozen atoms must be zero
    template <typename... Mems> void gradient(SimCell &, NeighbourList const &);

    // Compute gradient, force on frozen atoms must be zero
    template <typename... Mems> void hessian(SimCell const &, NeighbourList const &);

  private:
  };

}  // namespace otf