
// #pragma once

// #include <array>
// #include <cstddef>
// #include <optional>
// #include <vector>

// #include "libatom/asserts.hpp"
// #include "libatom/system/simbox.hpp"
// #include "libatom/utils.hpp"

// namespace otf {

//   /**
//    * @brief Maps {xyz} tuples to 1D grid.
//    *
//    * The Gridder is responisble for managing the mapping from nD->1D. Each atom is assigned to a
//    * cell (at least as large in each dimension as rcut) through .cell_idx(...). The neighbouring
//    * cells to any given cell can then be optained through a call to .neigh_cells(...).
//    */
//   class Gridder {
//   private:
//     floating m_rcut;
//     floating m_rcut_sq;

//     std::optional<OrthoSimBox> m_box;

//     Vec<int> m_shape;
//     Vec<int> m_prod_shape;

//     Vec<floating> m_cell;
//     Vec<floating> m_inv_cell;

//     std::vector<std::array<int, ipow<spatial_dims>(3) - 1>> m_neigh_cells;

//     /**
//      * @brief Cast a Vec<floating> to Vec<int>
//      */
//     Vec<int> to_ints(Vec<floating> const &x) const & {
//       ;
//       ASSERT(m_box, "Call before load");
//       return (x * m_inv_cell).cast<int>();
//     }

//     /**
//      * @brief Get the 1D index from the nD index
//      */
//     int to_1D(Vec<int> const &indexes) const & {
//       ;
//       ASSERT(m_box, "Call before load");
//       ASSERT((indexes > 0 && indexes < m_shape).all(), "Atom out of bounds");
//       return (indexes * m_prod_shape).sum();
//     }

//   public:
//     /**
//      * @brief Compute the cell index from an atoms position
//      *
//      * Assumes atom is in the cannonicle cell + displaced by m_cell
//      */
//     int cell_idx(Vec<floating> const &x) const {
//       ;

//       ASSERT((x >= 0).all(), "Atom out of bounds");

//       return to_1D(to_ints(x));
//     }

//     /**
//      * @brief Check each index, i, in each axis is in range: 1 <= i <= m_shape - 1
//      */
//     bool is_inner_cell(Vec<floating> const &x) {
//       ;

//       ASSERT(m_box, "Call before load");

//       Vec<int> indexes = to_ints(x);

//       return (indexes >= 1 && indexes <= m_shape - 1).all();
//     }

//     std::array<int, ipow<spatial_dims>(3) - 1> const &neigh_cells(std::size_t n) const noexcept {
//       return m_neigh_cells[n];
//     }

//     /**
//      * @brief Get the cell shape.
//      */
//     Vec<floating> const &origin() const noexcept { return m_cell; }

//     /**
//      * @brief Build the list of neighbour cells and memoize the results.
//      */
//     void compute_neigh_cells(OrthoSimBox const &box, double rcut);
//   };

// }  // namespace otf