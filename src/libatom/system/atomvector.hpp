#pragma once

#include <algorithm>
#include <array>
#include <cstddef>  // size_t
#include <cstdint>  // fast16_t
#include <optional>
#include <string_view>
#include <vector>

#include "aligned_allocator.hpp"
#include "bitsery/bitsery.h"
#include "bitsery/traits/array.h"
#include "bitsery/traits/vector.h"
#include "libatom/asserts.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief Represent a collection of atoms {x, z} stored as xxxx, zzzz.
   *
   */
  class AtomVector {
  private:
    static constexpr auto Align = Eigen::Aligned128;

  public:
    /**
     * @brief Get the number of atoms in the AtomVector.
     */
    [[nodiscard]] std::size_t size() const noexcept { return m_z.size(); }

    /**
     * @brief Fetch an uninitialised Mat3N large enough to hold the coordinates of every atom
     */
    [[nodiscard]] Mat3N<flt_t> empty_like_x() const {
      return {spatial_dims, static_cast<Eigen::Index>(size())};
    }

    /**
     * @brief Add an atom to the AtomVector
     *
     * @param x Coordinate of atom.
     * @param s Species number of the atom.
     */
    void emplace_back(Vec<flt_t> const& x, std::size_t z) {
      m_x.insert(m_x.end(), x.begin(), x.end());
      m_z.push_back(z);
    }

    /**
     * @brief Fetch an Eigen view into the atomic positions.
     */
    [[nodiscard]] Eigen::Map<Mat3N<flt_t>, Align> x() {
      STACK();
      ASSERT(m_x.size() % spatial_dims == 0, "Non integral number of atoms");
      return {m_x.data(), spatial_dims, static_cast<Eigen::Index>(m_x.size() / spatial_dims)};
    }

    /**
     * @brief Fetch a const Eigen view into the atomic positions.
     */
    [[nodiscard]] Eigen::Map<Mat3N<flt_t> const, Align> x() const {
      STACK();
      ASSERT(m_x.size() % spatial_dims == 0, "Non integral number of atoms");
      return {m_x.data(), spatial_dims, static_cast<Eigen::Index>(m_x.size() / spatial_dims)};
    }

    /**
     * @brief Fetch an Eigen view into the species numbers.
     */
    [[nodiscard]] Eigen::Map<VecN<std::size_t>, Align> z() {
      return {m_z.data(), static_cast<Eigen::Index>(m_z.size())};
    }

    /**
     * @brief Fetch a const Eigen view into the species numbers.
     */
    [[nodiscard]] Eigen::Map<VecN<std::size_t> const, Align> z() const {
      return {m_z.data(), static_cast<Eigen::Index>(m_z.size())};
    }

  private:
    std::vector<flt_t, detail::aligned<flt_t, Align>> m_x;
    std::vector<std::size_t, detail::aligned<std::size_t, Align>> m_z;

  protected:
    friend class bitsery::Access;

    /**
     * @brief Bitsery serialisation
     */
    template <typename S> void serialize(S& s) {
      STACK();

      ASSERT(size() < 1'000'000, "Too many atoms");

      s.container(m_x, 1'000'000 * 3);
      s.container(m_z, 1'000'000);
    }
  };

}  // namespace otf