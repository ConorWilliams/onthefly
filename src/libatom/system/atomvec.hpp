#pragma once

#include <algorithm>
#include <array>
#include <cstddef>  // size_t
#include <cstdint>  // fast16_t
#include <optional>
#include <string_view>
#include <vector>

#include "aligned_allocator.hpp"
#include "libatom/asserts.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief Represent atom type with <= 2-char Chemical symbol.
   */
  struct Symbol final : std::array<char, 2> {
    /**
     * @brief Construct a new Symbol object from a string
     */
    explicit Symbol(std::string_view sv) : std::array<char, 2>{} {
      STACK();
      ASSERT(sv.size() <= 2, "Symbols must be less than 2 chars");
      std::copy(sv.begin(), sv.end(), this->begin());
    }
  };

  /**
   * @brief Represent a collection of atoms {x, z} stored as xxxx, zzzz.
   */
  class AtomVector {
  private:
    static constexpr auto Align = Eigen::Aligned128;

    std::vector<Symbol> m_species_map;

    std::vector<flt_t, aligned<flt_t, Align>> m_x;
    std::vector<int_fast16_t, aligned<int_fast16_t, Align>> m_z;

  public:
    /**
     * @brief Get the number of atoms in the AtomVector.
     */
    std::size_t size() const noexcept { return m_z.size(); }

    /**
     * @brief Get the number of species in the AtomVector.
     */
    std::size_t num_species() const noexcept { return m_species_map.size(); }

    /**
     * @brief Convert a species value to the currently used integer to represent that species.
     *
     * @return std::optional<int_fast16_t> Empty if that species is not in the AtomVector.
     */
    std::optional<int_fast16_t> species2z(Symbol const& z) const noexcept {
      for (std::size_t i = 0; i < m_species_map.size(); i++) {
        if (z == m_species_map[i]) {
          return i;
        }
      }
      return {};
    }

    /**
     * @brief Add an atom to the AtomVector, If this is a new species updates the species map.
     *
     * @param x Coordinate of atom.
     * @param s Species of atom.
     */
    void emplace_back(Vec<flt_t> const& x, Symbol const& s) {
      // Insert position
      m_x.insert(m_x.end(), x.begin(), x.end());

      if (std::optional z = species2z(s)) {
        m_z.push_back(*z);
      } else {
        m_z.push_back(m_species_map.size());
        m_species_map.push_back(s);
      }
    }

    /**
     * @brief Fetch an Eigen view into the atomic positions.
     */
    Eigen::Map<Mat3N<flt_t>, Align> x() {
      STACK();
      ASSERT(m_x.size() % spatial_dims == 0, "Non integral number of atoms");
      return {m_x.data(), spatial_dims, static_cast<Eigen::Index>(m_x.size() / spatial_dims)};
    }

    /**
     * @brief Fetch a const Eigen view into the atomic positions.
     */
    Eigen::Map<Mat3N<flt_t> const, Align> x() const {
      STACK();
      ASSERT(m_x.size() % spatial_dims == 0, "Non integral number of atoms");
      return {m_x.data(), spatial_dims, static_cast<Eigen::Index>(m_x.size() / spatial_dims)};
    }

    /**
     * @brief Fetch an Eigen view into the species numbers.
     */
    Eigen::Map<VecN<int_fast16_t>, Align> z() {
      return {m_z.data(), static_cast<Eigen::Index>(m_z.size())};
    }

    /**
     * @brief Fetch a const Eigen view into the species numbers.
     */
    Eigen::Map<VecN<int_fast16_t> const, Align> z() const {
      return {m_z.data(), static_cast<Eigen::Index>(m_z.size())};
    }
  };

}  // namespace otf