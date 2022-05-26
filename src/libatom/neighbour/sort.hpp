#pragma once

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "libatom/atom_array.hpp"
#include "libatom/neighbour/grid.hpp"
#include "libatom/ortho_sim_box.hpp"

namespace otf::neighbour {

  /**
   * @brief Sort an atom array such that atoms in the same grid cell appear close to each other.
   *
   * @return AtomArray<Mems...> A sorted copy of the input.
   */
  template <typename... Mems> AtomArray<Mems...> sort(Grid grid, AtomArray<Mems...> const& in) {
    //
    struct Pair {
      std::size_t k;
      std::size_t j;
    };

    std::vector<Pair> aux(in.size());

    for (std::size_t i = 0; i < in.size(); i++) {
      aux[i].k = grid.cell_idx(grid.canon_grid_pos(in(Position{}, i)));
      aux[i].j = i;
    }

    std::sort(aux.begin(), aux.end(), [&](Pair const& a, Pair const& b) {
      return a.k == b.k ? in(Position{}, a.j)[0] < in(Position{}, b.j)[0] : a.k < b.k;
    });

    AtomArray<Mems...> out(in.size());

    for (std::size_t i = 0; i < in.size(); i++) {
      ((void)(out(Mems{}, i) = out(Mems{}, aux[i].j)), ...);
    }

    return out;
  }

}  // namespace otf::neighbour