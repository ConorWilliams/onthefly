#pragma once

#include <cstddef>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/atom_array.hpp"
#include "libatom/neighbour/grid.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::neighbour {

  /**
   * @brief The maximum nuber of ghosts neighbour cell supports is MAX_GHOST_RATIO * num_atoms
   */
  inline constexpr std::size_t MAX_GHOST_RATIO = 4;

  /**
   * @brief A class to contain, build and manage neighbour lists in shared memory.
   *
   * Designed with the intention of being reused List separates the building, updating and
   * using of neighbour lists (NLs). NLs can be constructed from a SimCell and then used to find all
   * atoms within some cut off of an atom efficiantly.
   *
   * List resolves periodicity using ghost atoms, these are stored and managed internally.
   *
   * An example of using a List to count the average number of atoms within rcut of each
   * atom:
   *
   * @code{.cpp}
   *
   * #include "libatom/neighbour/list.hpp"
   * #include "libatom/sim_cell.hpp"
   *
   * using namespace otf;
   *
   * SimCell atoms = ... // Initialise a set of atoms in a {10, 10, 10} cell.
   *
   * List nlist(atoms.box, 3.0);
   *
   * nlist.rebuild(atoms); // Build the NL.
   *
   * std::size_t num_neigh = 0;
   *
   * for (std::size_t i = 0; i < atoms.size(); i++) {
   *    nlist.for_neighbours(i, [&](std::size_t, floating, Vec3<floating> const&) {
   *        num_neigh++;
   *    });
   * }
   *
   * double avg_num_neigh = (double) num_neigh / atoms.size();
   *
   * @endcode
   */
  class List {
  public:
    /**
     * @brief Construct a new Neighbour List object. The cut off, rcut, must be smaller than the
     * minimum OrthSimCell extent.
     */
    List(OrthoSimBox const& box, floating rcut)
        : m_grid(box, rcut, true), m_rcut(rcut), m_rcut_sq(rcut * rcut) {}

    /**
     * @brief Build the internal neighbour lists in parallel with openMP.
     *
     * After a call to this function the for_neighbours methods can be used to iterate over all
     * atoms within rcut of any atom.
     */
    void rebuild(SimCell const& atoms, std::size_t num_threads);

    /**
     * @brief Update the positions of all atoms and ghosts to Position{} -= deltas
     *
     * Usefull if using a skin distance and wanting to avoid rebuilding the lists.
     */
    void update_positions(SimCell::underlying_t<Position> const& deltas);

    /**
     * @brief Call f(n, r, dr) for every neighbour n of atom i, within distance rcut.
     *
     * n is the neighbour index which could be a ghost or real atom, to convert to the index of the
     * real atom regardless use .image_to_real(n)
     *
     * dr is the minimum image vector joining i to n and r is the norm of dr
     */
    template <typename F> void for_neighbours(std::size_t i, floating rcut, F&& f) const {
      //
      ASSERT(rcut <= m_rcut, "Neighbour lists built with a larger rcut.");

      for (auto&& n : m_neigh_lists[i]) {
        Vec3<floating> const dr = m_atoms(Position{}, n) - m_atoms(Position{}, i);
        floating const r = norm(dr);
        if (r < rcut) {
          f(n, r, dr);
        }
      }
    }

    /**
     * @brief Call f(n, dr) for every neighbour n of atom i, within the cut off specified
     * during call to rebuild*.
     *
     * n is the neighbour index which could be a ghost or real atom, to convert to the index of the
     * real atom regardless use .image_to_real(n)
     *
     * dr is the minimum image vector joining i to n and r is the norm of dr
     */
    template <typename F> void for_neighbours(std::size_t i, F&& f) const {
      for (auto&& n : m_neigh_lists[i]) {
        f(n, m_atoms(Position{}, n) - m_atoms(Position{}, i));
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

    Grid m_grid;

    std::vector<std::size_t> m_head;

    floating m_rcut;

    floating m_rcut_sq;

    ///

    struct Next : AtomArrayMem<std::size_t, 1> {};

    AtomArray<Position, Index, Next> m_atoms;

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

}  // namespace otf::neighbour