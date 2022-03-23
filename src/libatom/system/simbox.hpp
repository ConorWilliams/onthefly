#pragma once

#include "libatom/asserts.hpp"
#include "libatom/utils.hpp"

namespace otf {
  /**
   * @brief Provides details of the simulations Orthoganal supercell geometry, all queries of the
   * space in which the atoms exist are provided by this class. It is assumed (and must be ensured)
   * all non-periodic atoms are within the OrthoSimBox extents.
   */
  class OrthoSimBox {
  public:
    /**
     * @brief Construct a new Ortho Sim Box object
     *
     * @param extents Length of simulation box along each axis.
     * @param periodic True for each periodic axis.
     */
    OrthoSimBox(Vec<flt_t> const &extents, Vec<bool> const &periodic);

    /**
     * @brief Maps atom into canonical cell, 0 <= r_i < extent_i for all i which are periodic.
     * Non-periodic atoms are within the simbox extents so v[i] * inv_extents less than 1 and  v[i]
     * remains unaffected, hence no non-periodic switch/select.
     */
    template <typename T>
    [[nodiscard]] Vec<flt_t> canon_image(Eigen::ArrayBase<T> const &v) const noexcept {
      STACK();
      ASSERT((m_periodic || (v >= Vec<flt_t>::Zero() && v < m_extents)).all(), "Out of box");
      return v - m_extents * (v * m_inv_extents).floor();
    }

    /**
     * @brief Compute the shortest vector connecting a to a periodic image of b. This function is
     * branchy and should be avoided in hot code.
     */
    template <typename A, typename B>
    [[nodiscard]] Vec<flt_t> min_image(Eigen::ArrayBase<A> const &a,
                                       Eigen::ArrayBase<B> const &b) const noexcept {
      Vec<flt_t> dr = b - a;
      return m_periodic.select(dr - m_extents * (dr * m_inv_extents + 0.5).floor(), dr);
    }

    friend bool operator==(OrthoSimBox const &a, OrthoSimBox const &b) noexcept {
      return (a.m_extents == b.m_extents).all() && (a.m_periodic == b.m_periodic).all();
    }

  private:
    Vec<flt_t> m_extents;
    Vec<bool> m_periodic;
    Vec<flt_t> m_inv_extents;
  };

}  // namespace otf