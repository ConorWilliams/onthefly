
#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/system/simbox.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief Maps {xyz} tuples to 1D grid.
   *
   * The Gridder is responisble for managing the mapping from nD to 1D. Each atom is assigned to a
   * cell (at least as large in each dimension as rcut) through .cell_idx(...). The neighbouring
   * cells to any given cell can then be optained through a call to .neigh_cells(...).
   */
  class Gridder {
  private:
    flt_t m_rcut;
    flt_t m_rcut_sq;

    std::optional<OrthoSimBox> m_box;

    Vec<int> m_shape;
    Vec<int> m_prod_shape;

    Vec<flt_t> m_cell;
    Vec<flt_t> m_inv_cell;

    std::vector<std::array<int, ipow<spatial_dims>(3) - 1>> m_neigh_cells;

  public:
    /**
     * @brief Compute the cell index from an atoms position
     *
     * Assumes atom is in the cannonicle cell + displaced by m_cell
     */
    int cell_idx(Vec<flt_t> const &x) const {
      STACK();

      ASSERT(m_box, "Call before load");

      Vec<int> indexes = (x * m_inv_cell).cast<int>();

      ASSERT((indexes >= 0).all(), "Atom out of bounds");
      ASSERT((indexes < m_shape).all(), "Atom out of bounds");

      return (indexes * m_prod_shape).sum();
    }

    int inner_cell_idx(Vec<flt_t> const &x) {
      STACK();

      ASSERT(m_box, "Call before load");

      Vec<int> indexes = (x * m_inv_cell).cast<int>();

      indexes
          = (indexes < 1)
                .select(Vec<int>::Ones(),
                        (indexes < m_shape - 1).select(indexes, (m_shape - 1) * Vec<int>::Ones()));

      return (indexes * m_prod_shape).sum();
    }

    std::array<int, ipow<spatial_dims>(3) - 1> const &neigh_cells(std::size_t n) const noexcept {
      return m_neigh_cells[n];
    }

    /**
     * @brief Get the cell shape.
     */
    Vec<flt_t> const &origin() const noexcept { return m_cell; }

    /**
     * @brief Build the list of neighbour cells and memoize the results.
     */
    void compute_neigh_cells(OrthoSimBox const &box, double rcut);
  };

}  // namespace otf