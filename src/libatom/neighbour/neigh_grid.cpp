#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/neighbour/neigh_grid.hpp"
#include "libatom/system/ortho_sim_box.hpp"
#include "libatom/utils.hpp"

namespace otf {

  NeighGrid::NeighGrid(OrthoSimBox const &box, floating rcut, bool compute_neigh_cells)
      : m_box{box}, m_rcut{rcut} {
    m_shape = 2 + (box.extents() / rcut).cast<int>();
    m_cell = box.extents() / (box.extents() / rcut).floor();
    m_inv_cell = 1.0 / m_cell;

    // Sanity checks
    VERIFY(rcut > 0, "rcut is negative");
    VERIFY((box.extents() > rcut).all(), "rcut is too big");

    // Cumprod _shape

    m_prod_shape = Vec3<int>::Ones();

    for (size_t i = 1; i < spatial_dims; i++) {
      m_prod_shape[i] = m_prod_shape[i - 1] * m_shape[i - 1];
    }

    if (compute_neigh_cells) {
      m_neigh_cells.resize(num_cells());

      for (int k = 0; k < m_shape[2]; k++) {
        for (int j = 0; j < m_shape[1]; j++) {
          for (int i = 0; i < m_shape[0]; i++) {
            //
            std::size_t next_slot = 1;

            // ^ Iterating over all {i, j, k} indexes in  the grid.
            for (int kk = std::max(0, k - 1); kk < std::min(k + 2, m_shape[2]); kk++) {
              for (int jj = std::max(0, j - 1); jj < std::min(j + 2, m_shape[1]); jj++) {
                for (int ii = std::max(0, i - 1); ii < std::min(i + 2, m_shape[0]); ii++) {
                  // ^ Iterating over every {ii, jj, kk} cell adjecent to {i, j, k}.
                  if (!(ii == i && jj == j && kk == k)) {
                    m_neigh_cells[to_1D({i, j, k})][next_slot++] = to_1D({ii, jj, kk});
                  }
                }
              }
            }
            // Write number of neigh_cells into first slot
            m_neigh_cells[to_1D({i, j, k})][0] = next_slot - 1;
          }
        }
      }
    }
  }

}  // namespace otf