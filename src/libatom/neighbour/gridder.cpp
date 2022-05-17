#include "libatom/neighbour/gridder.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/system/ortho_sim_box.hpp"
#include "libatom/utils.hpp"

namespace otf {

  void Gridder::compute_neigh_cells(OrthoSimBox const &box, double rcut) {
    // Early exit
    if (box == m_box && m_rcut == rcut) {
      return;
    }

    m_rcut = rcut;
    m_rcut_sq = rcut * rcut;

    m_box = box;

    m_shape = 2 + (box.extents() / rcut).cast<int>();
    m_cell = box.extents() / (box.extents() / rcut).floor();
    m_inv_cell = 1.0 / m_cell;

    // Sanity checks
    VERIFY(rcut > 0, "rcut is negative");
    VERIFY((box.extents() >= rcut).all(), "rcut is too big");

    // Cumprod _shape

    m_prod_shape = Vec3<int>::Ones();

    for (size_t i = 1; i < spatial_dims; i++) {
      m_prod_shape[i] = m_prod_shape[i - 1] * m_shape[i - 1];
    }

    std::array<int, ipow<spatial_dims>(3) - 1> offsets;

    // Compute neighbour stride offsets

    int idx = 0;

    for (auto k : {-1, 0, 1}) {
      for (auto j : {-1, 0, 1}) {
        if constexpr (spatial_dims == 3) {
          for (auto i : {-1, 0, 1}) {
            if (i != 0 || j != 0 || k != 0) {
              offsets[idx++] = (Vec3<int>{i, j, k} * m_prod_shape).sum();
            }
          }
        } else {
          if (j != 0 || k != 0) {
            offsets[idx++] = (Vec3<int>{j, k} * m_prod_shape).sum();
          }
        }
      }
    }

    // Fill in m_neigh_cells

    m_neigh_cells.resize(m_shape.prod());

    if constexpr (spatial_dims == 3) {
      for (int i = 1; i < m_shape[0] - 1; i++) {
        for (int j = 1; j < m_shape[1] - 1; j++) {
          for (int k = 1; k < m_shape[2] - 1; k++) {
            for (std::size_t n = 0; n < offsets.size(); n++) {
              auto lam = (Vec3<int>{i, j, k} * m_prod_shape).sum();
              m_neigh_cells[lam][n] = lam + offsets[n];
            }
          }
        }
      }
    }
  }

}  // namespace otf