#pragma once

#include "libatom/neighbour/neighbour_list.hpp"

namespace otf {

  void NeighbourCell::rebuild_neighbour_lists(SimCell const& atoms, double rcut) {
    // Attempt to avoid allocation
    if (atoms.size() < m_atoms.size()) {
      m_atoms.shrink(atoms.size());
    } else {
      m_atoms = AtomVector<Position, Index, Neighbours>(atoms.size());
    }

    m_num_atoms = atoms.size();

    m_grid.compute_neigh_cells(atoms.box, rcut);

    for (std::size_t i = 0; i < m_num_atoms; ++i) {
      _list.emplace_back(cell.canonicle_image(cell.activ[idx].vec) + _cell, cell.activ[idx].col);
    }

    make_ghosts();

    // Update head;
    _head.assign(_prod_shape[2], nullptr);

    for (auto& atom : _list) {
      atom.next = std::exchange(_head[lambda(atom)], &atom);
    }
  }

  void NeighbourCell::update_positions(SimCell const& atoms) {}

}  // namespace otf