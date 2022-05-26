
#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/ortho_sim_box.hpp"
#include "libatom/utils.hpp"
#include "nonstd/span.hpp"

namespace otf::neighbour {

  /**
   * @brief Maps {xyz} tuples to 1D grid.
   *
   * The Grid is responisble for managing the mapping from 3D->1D for neighbour lists. Each
   * atom is assigned to a cell (at least as large in each dimension as rcut) through
   * .cell_idx(...). The neighbouring cells to any given cell can then be optained through a call to
   * .neigh_cells(...).
   *
   * The canonicle cell is the cuboid of space spanned by the extents of the OrthoSimBox,
   * canon_grid_pos() maps a point into the canonicle cell and then nudges it by one grid cell in
   * the (1,1,1) direction. This leaves room for a layer of ghost atom cells around the canonicle
   * cell.
   */
  class Grid {
  public:
    /**
     * @brief Construct a new Neigh Grid object, optionally chose to build the internal list of
     * neighbour cells required for calls to neigh_cells().
     */
    Grid(OrthoSimBox const &box, floating rcut, bool compute_neigh_cells);

    /**
     * @brief Get shape of a single cell
     */
    Vec3<floating> const &cell() const noexcept { return m_cell; }

    /**
     * @brief Get the total number of cells in the neighbour grid.
     */
    std::size_t num_cells() const noexcept { return m_shape.prod(); }

    /**
     * @brief Compute the cell index from an atoms position
     *
     * Assumes atom is in the cannonicle cell + displaced by m_cell
     */
    std::size_t cell_idx(Vec3<floating> const &x) const { return to_1D(clamp_to_grid_idxs(x)); }

    /**
     * @brief Get the Canonicle grid position, this is the canonicle image displaced by one grid
     * cell.
     */
    template <typename E>
    Vec3<floating> canon_grid_pos(Eigen::ArrayBase<E> const &x) const noexcept {
      return m_box.canon_image(x) + m_cell;
    }

    /**
     * @brief Get an span containing the indexes of every cell adjecent to the nth cell.
     *
     * Does not include the nth cell.
     */
    nonstd::span<std::size_t const> neigh_cells(std::size_t n) const {
      ASSERT(n < m_neigh_cells.size(), "Forgot to calc neigh cells");
      return {m_neigh_cells[n].data() + 1, m_neigh_cells[n][0]};
    }

  private:
    // Always ready

    Vec3<int> m_shape = Vec3<int>::Zero();
    Vec3<int> m_prod_shape = Vec3<int>::Zero();

    Vec3<floating> m_cell = Vec3<floating>::Zero();
    Vec3<floating> m_inv_cell = Vec3<floating>::Zero();

    // Memoization

    // Maximum of 3^3 - 1 neighbour cells, use first slot to store count of neighbours.
    std::vector<std::array<std::size_t, ipow<spatial_dims>(3)>> m_neigh_cells;

    OrthoSimBox m_box;

    floating m_rcut = 0;

    /**
     * @brief Cast a Vec3<floating> to Vec3<int> and clammp on interval [0, m_shape -1]
     */
    Vec3<int> clamp_to_grid_idxs(Vec3<floating> const &x) const & {
      return (x * m_inv_cell).cast<int>().cwiseMax(0).cwiseMin(m_shape - 1);
    }

    /**
     * @brief Get the 1D index from the nD index
     */
    std::size_t to_1D(Vec3<int> const &indexes) const & { return (indexes * m_prod_shape).sum(); }
  };

}  // namespace otf::neighbour