#pragma once

#include <cmath>

#include "libatom/asserts.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief Provides details of the simulations Orthoganal supercell geometry,
   *
   * All queries of the space in which the atoms exist are provided by this class. It is assumed
   * (and must be ensured) all non-periodic atoms are within the OrthoSimBox extents.
   */
  class OrthoSimBox {
  public:
    /**
     * @brief Construct an empty Ortho Sim Box object.
     */
    OrthoSimBox() = default;

    /**
     * @brief Construct a new Ortho Sim Box object.
     *
     * @param extents Length of simulation box along each axis.
     * @param periodic True for each periodic axis.
     */
    OrthoSimBox(Vec3<floating> const &extents, Vec3<bool> const &periodic)
        : m_extents{extents}, m_periodic{periodic}, m_inv_extents(1.0 / extents) {
      VERIFY((m_extents > 0).all(), "OrthoSimBox extents are negative");
    }

    /**
     * @brief Extents getter (const).
     */
    [[nodiscard]] Vec3<floating> const &extents() const noexcept { return m_extents; }

    /**
     * @brief Maps atom into canonical cell, 0 <= r_i < extent_i for all i which are periodic.
     * Non-periodic atoms are within the simbox extents so x[i] * inv_extents less than 1 and x[i]
     * remains unaffected, hence no non-periodic switch/select.
     */
    template <typename T>
    [[nodiscard]] Vec3<floating> canon_image(Eigen::ArrayBase<T> const &x) const noexcept {
      ASSERT((m_periodic || (x >= Vec3<floating>::Zero() && x < m_extents)).all(), "Out of box");
      return x - m_extents * (x * m_inv_extents).floor();
    }

    /**
     * @brief Compute the shortest Vector connecting a to a periodic image of b. This function is
     * branchy and should be avoided in hot code.
     */
    template <typename A, typename B>
    [[nodiscard]] Vec3<floating> min_image(Eigen::ArrayBase<A> const &a,
                                           Eigen::ArrayBase<B> const &b) const noexcept {
      Vec3<floating> dr = b - a;
      return m_periodic.select(dr - m_extents * (dr * m_inv_extents + 0.5).floor(), dr);
    }

    [[nodiscard]] friend bool operator==(OrthoSimBox const &a, OrthoSimBox const &b) noexcept {
      return (a.m_extents == b.m_extents).all() && (a.m_periodic == b.m_periodic).all();
    }

  private:
    Vec3<floating> m_extents;
    Vec3<bool> m_periodic;
    Vec3<floating> m_inv_extents;
  };

}  // namespace otf