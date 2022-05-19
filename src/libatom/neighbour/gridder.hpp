
#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/system/ortho_sim_box.hpp"
#include "libatom/utils.hpp"
#include "nonstd/span.hpp"

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
    int cell_idx(Vec3<floating> const &x) const { return to_1D(clamp_to_grid_idxs(x)); }

    /**
     * @brief Get an span containing the indexes of every cell adjecent to the nth cell.
     *
     * Does not include the nth cell.
     */
    nonstd::span<std::size_t const> neigh_cells(std::size_t n) const noexcept {
      return {m_neigh_cells[n].data() + 1, m_neigh_cells[n][0]};
    }

    /**
     * @brief Get the Canonicle grid position, this is the canonicle image displaced by one grid
     * cell.
     */
    template <typename E>
    Vec3<floating> canon_grid_pos(Eigen::ArrayBase<E> const &x) const noexcept {
      return m_box.canon_image(x) + m_cell;
    }

    /**
     * @brief Get shape of a single cell
     */
    Vec3<floating> const &cell() const noexcept { return m_cell; }

    std::size_t num_cells() const noexcept { return m_shape.prod(); }

    /**
     * @brief Build the list of neighbour cells and memoize the results.
     */
    void compute_neigh_cells(OrthoSimBox const &box, floating rcut);

  private:
    floating m_rcut = 0;
    floating m_rcut_sq = 0;

    OrthoSimBox m_box;

    Vec3<int> m_shape = Vec3<int>::Zero();
    Vec3<int> m_prod_shape = Vec3<int>::Zero();

    Vec3<floating> m_cell = Vec3<floating>::Zero();
    Vec3<floating> m_inv_cell = Vec3<floating>::Zero();

    // Maximum of 3^3 - 1 neighbour cells, use first slot to store count of neighbours.
    std::vector<std::array<std::size_t, ipow<spatial_dims>(3)>> m_neigh_cells;

    /**
     * @brief Cast a Vec3<floating> to Vec3<int> and clammp on interval [0, m_shape -1]
     */
    Vec3<int> clamp_to_grid_idxs(Vec3<floating> const &x) const & {
      return (x * m_inv_cell).cast<int>().cwiseMax(0).cwiseMin(m_shape - 1);
    }

    /**
     * @brief Get the 1D index from the nD index
     */
    int to_1D(Vec3<int> const &indexes) const & { return (indexes * m_prod_shape).sum(); }
  };

}  // namespace otf