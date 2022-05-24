#pragma once

#include <cstddef>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/neighbour/neigh_grid.hpp"
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
   * @brief A class to contain, build and manage neighbour lists in shared memory.
   *
   * Designed with the intention of being reused NeighbourList separates the building, updating and
   * using of neighbour lists (NLs). NLs can be constructed from a SimCell and then used to find all
   * atoms within some cut off of an atom efficiantly.
   *
   * NeighbourList resolves periodicity using ghost atoms, these are stored and managed internally.
   *
   * An example of using a NeighbourList to count the average number of atoms within rcut of each
   * atom:
   *
   * @code{.cpp}
   *
   * #include "libatom/neighbour/neighbour_list.hpp"
   * #include "libatom/system/sim_cell.hpp"
   *
   * using namespace otf;
   *
   * SimCell atoms = ... // Initialise a set of atoms in a {10, 10, 10} cell.
   *
   * NeighbourList nlist(atoms.box, 3.0);
   *
   * nlist.rebuild(atoms); // Build the NL.
   *
   * std::size_t num_neigh = 0;
   *
   * for (size_t i = 0; i < atoms.size(); i++) {
   *    nlist.for_neighbours(i, [&](std::size_t, floating, Vec3<floating> const&) {
   *        num_neigh++;
   *    });
   * }
   *
   * double avg_num_neigh = (double) num_neigh / atoms.size();
   *
   * @endcode
   */
  class NeighbourList {
  public:
    /**
     * @brief Construct a new Neighbour List object. The cut off, rcut, must be smaller than the
     * minimum OrthSimCell extent.
     */
    NeighbourList(OrthoSimBox const& box, floating rcut)
        : m_grid(box, rcut, true), m_rcut(rcut), m_rcut_sq(rcut * rcut) {}

    /**
     * @brief Build the internal neighbour lists in parallel with openMP.
     *
     * After a call to this function the for_neighbours methods can be used to iterate over all
     * atoms within rcut of any atom.
     */
    void rebuild(SimCell const& atoms, std::size_t num_threads = 1);

    /**
     * @brief Update the positions of all atoms + ghosts but do not rebuild the neighbour_lists.
     *
     * Usefull if using a skin distance.
     */
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
     * during call to rebuild*.
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
    //

    NeighGrid m_grid;

    std::vector<std::size_t> m_head;

    floating m_rcut;

    floating m_rcut_sq;

    ///

    struct Next : Member<std::size_t, 1> {};
    struct Offset : Member<floating, 3> {};

    AtomArray<Position, Index, Next, Offset> m_atoms;

    std::size_t m_num_plus_ghosts = 0;

    std::vector<std::vector<std::size_t>> m_neigh_lists;

    /**
     * @brief Initialise memory, load in atom atoms, build ghosts and kint link cell lists.
     */
    void init_and_build_lcl(SimCell const& atoms);

    /**
     * @brief Set up all ghost indexes positions and offsets.
     *
     * A ghosts position can be calculated from the position of its image plus its offset.
     */
    void make_ghosts(OrthoSimBox const& box);

    /**
     * @brief Build the neighbour list of the ith atom.
     */
    void build_neigh_list(std::size_t i);
  };

}  // namespace otf