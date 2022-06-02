#pragma once

#include <cmath>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
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
     * @brief Periodicity getter (const).
     */
    [[nodiscard]] Vec3<bool> const &periodic() const noexcept { return m_periodic; }

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
    Vec3<floating> m_extents = Vec3<floating>::Zero();
    Vec3<bool> m_periodic = Vec3<bool>::Zero();
    Vec3<floating> m_inv_extents = Vec3<floating>::Zero();
  };

  /**
   * @brief A SimCell is a collection of atoms augmented with an OrthoSimBox
   */
  class SimCell : public OrthoSimBox,
                  public AtomArray<Frozen, AtomicNum, Gradient, Position, Axis> {
  public:
    SimCell(OrthoSimBox const &arg) : OrthoSimBox{arg} {}

    void remove_soft_modes() {
      if (count_frozen() == 0) {
        for (std::size_t i = 0; i < spatial_dims; i++) {
          if (this->periodic()[i]) {
            (*this)(Gradient{}).row(i) = remove_soft_mode((*this)(Gradient{}).row(i));
          }
        }
      }
    }

    /**
     * @brief Count the number of frozen atoms in the SimCell
     */
    [[nodiscard]] std::size_t count_frozen() const { return (*this)(Frozen{}).count(); }

    /**
     * @brief Allocate (if required) and zero the storage for the hessian matrix.
     *
     * The hessian is an 3m by 3m matrix with m the number of active atoms. Each time this function
     * is called it will check the number of active atoms has not changed, if it has it will
     * reallocate a new hessian matrix *destroying* the previous one in the process.
     */
    void zero_hess() {
      std::size_t m = size() - count_frozen();
      m_hess = Eigen::Array<floating, Eigen::Dynamic, Eigen::Dynamic>::Zero(3 * m, 3 * m);
    }

    /**
     * @brief Get the full hessian matrix.
     */
    [[nodiscard]] auto hess() { return m_hess.topLeftCorner(m_hess.cols(), m_hess.rows()); }

    /**
     * @brief Get the spatial_dims * spatial_dims block of the hessian corresponding to the a^th and
     * b^th atom.
     */
    [[nodiscard]] auto hess(std::size_t a, std::size_t b) {
      return m_hess.block<spatial_dims, spatial_dims>(spatial_dims * a, spatial_dims * b);
    }

  private:
    Eigen::Array<floating, Eigen::Dynamic, Eigen::Dynamic> m_hess;
  };

}  // namespace otf