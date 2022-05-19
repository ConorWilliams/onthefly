#pragma once

#include <cstddef>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/neighbour/gridder.hpp"
#include "libatom/system/atom_array.hpp"
#include "libatom/system/member.hpp"
#include "libatom/system/ortho_sim_box.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief The maximum nuber of ghosts neighbour cell supports is MAX_GHOST_RATIO * num_atoms
   */
  inline constexpr std::size_t MAX_GHOST_RATIO = 4;
  /**
   * @brief
   *
   * Designed with the intention of reusing
   *
   */
  class NeighbourCell {
  public:
    void rebuild_neighbour_lists(SimCell const& atoms, floating rcut);

    void rebuild_neighbour_lists_parallel(SimCell const& atoms, floating rcut);

    void update_positions(SimCell const& atoms);

    /**
     * @brief Call f(n, r_sq, dr) for every neighbour n of atom i, within distance rcut.
     *
     * n is the neighbour index which could be a ghost or real atom, to convert to the index of the
     * real atom regardless of weather it is a ghost or not call .image_to_real(n)
     *
     * dr is the minimum image vector joining i to n and r_sq is the norm_sq of dr
     */
    template <typename F> void for_neighbours(std::size_t i, floating rcut, F&& f) const {
      //
      ASSERT(rcut < m_rcut, "Neighbour lists built with a larger rcut.");

      for (auto&& n : m_neigh_lists[i]) {
        Vec3<floating> const dr = m_atoms(Position{}, n) - m_atoms(Position{}, i);
        floating const r_sq = norm_sq(dr);
        if (r_sq < rcut * rcut) {
          f(n, r_sq, dr);
        }
      }
    }

    /**
     * @brief Call f(n, r_sq, dr) for every neighbour n of atom i, within the cut off specified
     * during call to rebuild_neighbour_lists*.
     *
     * n is the neighbour index which could be a ghost or real atom, to convert to the index of the
     * real atom regardless of weather it is a ghost or not call .image_to_real(n)
     *
     * dr is the minimum image vector joining i to n and r_sq is the norm_sq of dr
     */
    template <typename F> void for_neighbours(std::size_t i, F&& f) const {
      for (auto&& n : m_neigh_lists[i]) {
        Vec3<floating> const dr = m_atoms(Position{}, n) - m_atoms(Position{}, i);
        f(n, norm_sq(dr), dr);
      }
    }

    /**
     * @brief Convert the neighbour index of a real or ghost atom to the index of the real atom.
     *
     * This is separate, rather than the n being provided during for_neighbours so as not to pay for
     * a cache miss if the real index is not required.
     */
    std::size_t image_to_real(std::size_t i) const { return m_atoms(Index{}, i); }

  private:
    struct Next : Member<std::size_t, 1> {};

    struct Offset : Member<floating, 3> {};

    Gridder m_grid;

    AtomArray<Position, Index, Next, Offset> m_atoms;

    std::size_t m_num_plus_ghosts = 0;

    std::vector<std::vector<std::size_t>> m_neigh_lists;

    std::vector<std::size_t> m_head;

    floating m_rcut = 0;

    /**
     * @brief Initialise memory, load in atom atoms, build ghosts and kint link cell lists.
     */
    void init_and_build_lcl(SimCell const& atoms, floating rcut);

    /**
     * @brief Set up all ghost indexes positions and offsets.
     *
     * A ghosts position can be calculated from the position of its image plus its offset.
     */
    void make_ghosts(OrthoSimBox const& box, floating rcut);

    /**
     * @brief Broadcast new positions to ghosts.
     */
    void update_ghosts();

    /**
     * @brief Build the neighbour list of the ith atom.
     */
    void build_neigh_list(std::size_t i, floating rcut);
  };

}  // namespace otf