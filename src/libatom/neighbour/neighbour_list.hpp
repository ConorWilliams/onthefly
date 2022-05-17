#pragma once

#include <cstddef>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/neighbour/gridder.hpp"
#include "libatom/system/atom_array.hpp"
#include "libatom/system/member.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief The maximum nuber of ghosts neighbour cell supports is MAX_GHOST_RATIO * num_atoms
   */
  inline constexpr std::size_t MAX_GHOST_RATIO = 9;

  class NeighbourCell {
  public:
    void rebuild_neighbour_lists(SimCell const& atoms, double rcut);

    void update_positions(SimCell const& atoms);

  private:
    Gridder m_grid;

    AtomArray<Position, Index> m_atoms;

    std::vector<std::vector<std::size_t>> m_neigh_lists;
  };

}  // namespace otf