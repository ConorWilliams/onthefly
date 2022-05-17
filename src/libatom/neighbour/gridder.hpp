
#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/system/ortho_sim_box.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief Maps {xyz} tuples to 1D grid.
   *
   * The Gridder is responisble for managing the mapping from 3D->1D. Each atom is assigned to a
   * cell (at least as large in each dimension as rcut) through .cell_idx(...). The neighbouring
   * cells to any given cell can then be optained through a call to .neigh_cells(...).
   *
   * The canonicle cell is the cuboid of space spanned by the extents of the OrthoSimBox
   */
  class Gridder {
  public:
    /**
     * @brief Compute the cell index from an atoms position
     *
     * Assumes atom is in the cannonicle cell + displaced by m_cell
     */
    int cell_idx(Vec3<floating> const &x) const {
      ASSERT((x >= 0).all(), "Atom out of bounds");
      return to_1D(to_ints(x));
    }

    /**
     * @brief Check each index, i, in each axis is in range: 1 <= i <= m_shape - 1
     */
    bool is_inner_cell(Vec3<floating> const &x) {
      Vec3<int> indexes = to_ints(x);
      return (indexes >= 1 && indexes <= m_shape - 1).all();
    }

    std::array<int, ipow<spatial_dims>(3) - 1> const &neigh_cells(std::size_t n) const noexcept {
      return m_neigh_cells[n];
    }

    /**
     * @brief Get the cell shape.
     */
    Vec3<floating> const &cell() const noexcept { return m_cell; }

    /**
     * @brief Build the list of neighbour cells and memoize the results.
     */
    void compute_neigh_cells(OrthoSimBox const &box, double rcut);

  private:
    floating m_rcut;
    floating m_rcut_sq;

    OrthoSimBox m_box;

    Vec3<int> m_shape;
    Vec3<int> m_prod_shape;

    Vec3<floating> m_cell;
    Vec3<floating> m_inv_cell;

    std::vector<std::array<int, ipow<spatial_dims>(3) - 1>> m_neigh_cells;

    /**
     * @brief Cast a Vec3<floating> to Vec3<int>
     */
    Vec3<int> to_ints(Vec3<floating> const &x) const & { return (x * m_inv_cell).cast<int>(); }

    /**
     * @brief Get the 1D index from the nD index
     */
    int to_1D(Vec3<int> const &indexes) const & {
      ASSERT((indexes > 0 && indexes < m_shape).all(), "Atom out of bounds");
      return (indexes * m_prod_shape).sum();
    }
  };

}  // namespace otf