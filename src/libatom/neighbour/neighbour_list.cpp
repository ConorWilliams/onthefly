

#include "libatom/neighbour/neighbour_list.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <utility>

#include "libatom/asserts.hpp"
#include "libatom/system/member.hpp"
#include "libatom/system/ortho_sim_box.hpp"
#include "libatom/utils.hpp"

namespace otf {

  void NeighbourCell::rebuild_neighbour_lists(SimCell const& atoms, floating rcut) {
    //
    m_rcut = rcut;

    //
    if (m_neigh_lists.size() != atoms.size()) {
      // Allocate space if needed
      m_atoms.destructive_resize(atoms.size() * (1 + MAX_GHOST_RATIO));
      m_neigh_lists.resize(atoms.size(), {});
    }

    m_grid.compute_neigh_cells(atoms.box, rcut);

    // Copy in atoms

    for (std::size_t i = 0; i < atoms.size(); ++i) {
      m_atoms(Position{}, i) = atoms.box.canon_image(atoms(Position{}, i)) + m_grid.cell();
      m_atoms(Index{}, i) = i;
    }

    make_ghosts(atoms.box, rcut);

    // Update head.
    m_head.assign(m_grid.num_cells(), std::numeric_limits<std::size_t>::max());

    // Build LCL.
    for (size_t i = 0; i < m_num_plus_ghosts; i++) {
      m_atoms(Next{}, i) = std::exchange(m_head[m_grid.cell_idx(m_atoms(Position{}, i))], i);
    }

    build_neigh_lists(rcut);
  }

  void NeighbourCell::update_positions(SimCell const&) {}

  void NeighbourCell::build_neigh_lists(floating rcut) {
    floating rcut_sq = rcut * rcut;

    for (size_t i = 0; i < m_neigh_lists.size(); i++) {
      //
      m_neigh_lists[i].clear();

      std::size_t i_cell = m_grid.cell_idx(m_atoms(Position{}, i));

      std::size_t n = m_head[i_cell];

      // In same cell must check not-self
      while (n != std::numeric_limits<std::size_t>::max()) {
        if (n != i) {
          if (norm_sq(m_atoms(Position{}, i) - m_atoms(Position{}, n)) < rcut_sq) {
            m_neigh_lists[i].push_back(n);
          }
        }
        n = m_atoms(Next{}, n);
      }

      // In adjacent cells -- don't need to check against self
      for (auto&& n_cell : m_grid.neigh_cells(i_cell)) {
        //
        std::size_t n = m_head[n_cell];

        while (n != std::numeric_limits<std::size_t>::max()) {
          if (norm_sq(m_atoms(Position{}, i) - m_atoms(Position{}, n)) < rcut_sq) {
            m_neigh_lists[i].push_back(n);
          }
          n = m_atoms(Next{}, n);
        }
      }
    }
  }

  void NeighbourCell::make_ghosts(OrthoSimBox const& box, floating rcut) {
    //
    std::size_t next_slot = m_neigh_lists.size();

    for (std::size_t i = 0; i < spatial_dims; ++i) {
      // Only make ghosts if axis is periodic
      if (box.periodic()[i]) {
        //
        std::size_t const end = next_slot;

        for (std::size_t j = 0; j < end; ++j) {
          if (m_atoms(Position{}, j)[i] < m_grid.cell()[i] + rcut) {
            //
            std::size_t slot = next_slot++;

            ASSERT(slot < m_atoms.size(), "Not enough space for ghosts");

            m_atoms(Position{}, slot) = m_atoms(Position{}, j);
            m_atoms(Position{}, slot)[i] += box.extents()[i];
            m_atoms(Index{}, slot) = m_atoms(Index{}, j);
          }

          if (m_atoms(Position{}, j)[i] >= box.extents()[i] + m_grid.cell()[i] - rcut) {
            //
            std::size_t slot = next_slot++;

            ASSERT(slot < m_atoms.size(), "Not enough space for ghosts");

            m_atoms(Position{}, slot) = m_atoms(Position{}, j);
            m_atoms(Position{}, slot)[i] -= box.extents()[i];
            m_atoms(Index{}, slot) = m_atoms(Index{}, j);
          }
        }
      }
    }

    m_num_plus_ghosts = next_slot;
  }

}  // namespace otf