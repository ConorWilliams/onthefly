

#include "libatom/neighbour/neighbour_list.hpp"

#include "libatom/system/member.hpp"

namespace otf {

  void NeighbourCell::rebuild_neighbour_lists(SimCell const& atoms, double rcut) {
    //
    if (m_neigh_lists.size() != atoms.size()) {
      // Allocate space if needed
      m_atoms.destructive_resize(atoms.size() * (1 + MAX_GHOST_RATIO));
      m_neigh_lists.resize(atoms.size(), {});
    }

    m_grid.compute_neigh_cells(atoms.box, rcut);

    for (std::size_t i = 0; i < atoms.size(); ++i) {
      m_atoms(Position{}, i) = atoms.box.canon_image(atoms(Position{}, i)) + m_grid.cell();

      //   _list.emplace_back(, cell.activ[idx].col);
    }

    // make_ghosts();

    // // Update head;
    // _head.assign(_prod_shape[2], nullptr);

    // for (auto& atom : _list) {
    //   atom.next = std::exchange(_head[lambda(atom)], &atom);
    // }
  }

  void NeighbourCell::update_positions(SimCell const&) {}

}  // namespace otf