// #pragma once

// #include <cmath>

// #include "bitsery/bitsery.h"
// #include "libatom/asserts.hpp"
// #include "libatom/utils.hpp"

// namespace otf {
//   /**
//    * @brief Provides details of the simulations Orthoganal supercell geometry,
//    *
//    * All queries of the space in which the atoms exist are provided by this class. It is assumed
//    * (and must be ensured) all non-periodic atoms are within the OrthoSimBox extents.
//    */
//   class OrthoSimBox {
//   public:
//     /**
//      * @brief Construct a new Ortho Sim Box object
//      *
//      * @param extents Length of simulation box along each axis.
//      * @param periodic True for each periodic axis.
//      */
//     OrthoSimBox(Vec<floating> const &extents, Vec<bool> const &periodic);

//     /**
//      * @brief Extents getter (const)
//      */
//     Vec<floating> const &extents() const noexcept { return m_extents; }

//     /**
//      * @brief Maps atom into canonical cell, 0 <= r_i < extent_i for all i which are periodic.
//      * Non-periodic atoms are within the simbox extents so x[i] * inv_extents less than 1 and
//      x[i]
//      * remains unaffected, hence no non-periodic switch/select.
//      */
//     template <typename T>
//     [[nodiscard]] Vec<floating> canon_image(Eigen::ArrayBase<T> const &x) const noexcept {
//       ;
//       ASSERT((m_periodic || (x >= Vec<floating>::Zero() && x < m_extents)).all(), "Out of box");
//       return x - m_extents * (x * m_inv_extents).floor();
//     }

//     /**
//      * @brief Compute the shortest vector connecting a to a periodic image of b. This function is
//      * branchy and should be avoided in hot code.
//      */
//     template <typename A, typename B>
//     [[nodiscard]] Vec<floating> min_image(Eigen::ArrayBase<A> const &a,
//                                           Eigen::ArrayBase<B> const &b) const noexcept {
//       Vec<floating> dr = b - a;
//       return m_periodic.select(dr - m_extents * (dr * m_inv_extents + 0.5).floor(), dr);
//     }

//     friend bool operator==(OrthoSimBox const &a, OrthoSimBox const &b) noexcept {
//       return (a.m_extents == b.m_extents).all() && (a.m_periodic == b.m_periodic).all();
//     }

//     /**
//      * @brief Not part of the API, for internal use only
//      *
//      * @return OrthoSimBox Unitialised!
//      */
//     [[nodiscard]] static OrthoSimBox empty() { return {}; }

//   private:
//     Vec<floating> m_extents;
//     Vec<bool> m_periodic;
//     Vec<floating> m_inv_extents;

//   protected:
//     friend class bitsery::Access;

//     /**
//      * @brief Construct a new OrthoSimBox object don't worry about class invariants, they will be
//      * restored in deserialization
//      */
//     OrthoSimBox() = default;

//     /**
//      * @brief Bitsery serialisation
//      */
//     template <typename S> void serialize(S &s) { s(m_extents, m_periodic, m_inv_extents); }
//   };

// }  // namespace otf