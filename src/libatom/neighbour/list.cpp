

#include "libatom/neighbour/list.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>

#include "libatom/asserts.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::neighbour {

  void List::rebuild(SimCell const& atoms, std::size_t num_threads) {
    //
    init_and_build_lcl(atoms);

#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (std::size_t i = 0; i < m_neigh_lists.size(); i++) {
      build_neigh_list(i);
    }
  }

  void List::init_and_build_lcl(SimCell const& atoms) {
    //
    if (m_neigh_lists.size() != atoms.size()) {
      // Allocate space if needed
      m_atoms.destructive_resize(atoms.size() * (1 + MAX_GHOST_RATIO));
      m_neigh_lists.resize(atoms.size(), {});

      // Need to re-index if size changed
      for (std::size_t i = 0; i < atoms.size(); ++i) {
        m_atoms(Index{}, i) = i;
      }
    }

    // Copy in atoms
    for (std::size_t i = 0; i < atoms.size(); ++i) {
      m_atoms(Position{}, i) = m_grid.canon_grid_pos(atoms(Position{}, i));
    }

    make_ghosts(atoms.box);

    // Update head.
    m_head.assign(m_grid.num_cells(), std::numeric_limits<std::size_t>::max());

    // Build LCL.
    for (std::size_t i = 0; i < m_num_plus_ghosts; i++) {
      m_atoms(Next{}, i) = std::exchange(m_head[m_grid.cell_idx(m_atoms(Position{}, i))], i);
    }
  }

  void List::update_positions(SimCell::underlying_t<Position> const& deltas) {
    // Copy in atoms
    for (std::size_t i = 0; i < m_neigh_lists.size(); ++i) {
      m_atoms(Position{}, i) -= deltas.col(i);
    }

    // Update ghosts
    for (std::size_t i = m_neigh_lists.size(); i < m_num_plus_ghosts; i++) {
      m_atoms(Position{}, i) -= deltas.col(image_to_real(i));
    }
  }

  void List::build_neigh_list(std::size_t i) {
    //
    m_neigh_lists[i].clear();

    std::size_t i_cell = m_grid.cell_idx(m_atoms(Position{}, i));

    {
      std::size_t n = m_head[i_cell];

      // In same cell must check not-self
      while (n != std::numeric_limits<std::size_t>::max()) {
        if (n != i) {
          if (norm(m_atoms(Position{}, i) - m_atoms(Position{}, n)) < m_rcut) {
            m_neigh_lists[i].push_back(n);
          }
        }
        n = m_atoms(Next{}, n);
      }
    }

    // In adjacent cells -- don't need to check against self
    for (auto&& n_cell : m_grid.neigh_cells(i_cell)) {
      //
      std::size_t n = m_head[n_cell];

      while (n != std::numeric_limits<std::size_t>::max()) {
        if (norm(m_atoms(Position{}, i) - m_atoms(Position{}, n)) < m_rcut) {
          m_neigh_lists[i].push_back(n);
        }
        n = m_atoms(Next{}, n);
      }
    }
  }

  void List::make_ghosts(OrthoSimBox const& box) {
    //
    std::size_t next_slot = m_neigh_lists.size();

    for (std::size_t i = 0; i < spatial_dims; ++i) {
      // Only make ghosts if axis is periodic
      if (box.periodic()[i]) {
        //
        std::size_t const end = next_slot;

        for (std::size_t j = 0; j < end; ++j) {
          if (m_atoms(Position{}, j)[i] < m_grid.cell()[i] + m_rcut) {
            //
            std::size_t slot = next_slot++;

            ASSERT(slot < m_atoms.size(), "Not enough space for ghosts");

            m_atoms(Position{}, slot) = m_atoms(Position{}, j);
            m_atoms(Position{}, slot)[i] += box.extents()[i];

            m_atoms(Index{}, slot) = m_atoms(Index{}, j);
          }

          if (m_atoms(Position{}, j)[i] >= box.extents()[i] + m_grid.cell()[i] - m_rcut) {
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

}  // namespace otf::neighbour