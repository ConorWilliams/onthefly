#pragma once

#include <cstddef>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/neighbour/gridder.hpp"
#include "libatom/system/atom_vector.hpp"
#include "libatom/system/member.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf {

  struct Neighbours : Member<std::vector<std::size_t>, 1> {};

  class NeighbourCell {
  public:
    void rebuild_neighbour_lists(SimCell const& atoms, double rcut);

    void update_positions(SimCell const& atoms);

  private:
    Gridder m_grid;

    AtomVector<Position, Index, Neighbours> m_atoms;

    std::size_t m_num_atoms = 0;
  };

}  // namespace otf