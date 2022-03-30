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
   * @brief Represent atom type with a 1/2-char Chemical symbol.
   */
  struct Symbol final : std::array<char, 2> {
    /**
     * @brief Construct a new Symbol object from a string, pad with char{}.
     */
    explicit Symbol(std::string_view sv) : std::array<char, 2>{} {
      STACK();
      ASSERT(sv.size() == 1 || sv.size() == 2, "Symbols must be less than 2 chars");
      std::copy(sv.begin(), sv.end(), this->begin());
    }

  private:
    friend class bitsery::Access;

    /**
     * @brief Construct a new OrthoSimBox object don't worry about class invariants, they will be
     * restored in deserialization
     */
    Symbol() = default;

    /**
     * @brief Bitsery serialisation
     */
    template <typename S> void serialize(S& s) {
      s.container(static_cast<std::array<char, 2>&>(*this));
    }
  };

  /**
   * @brief Represent a collection of atoms {x, z} stored as xxxx, zzzz.
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
     * @brief Get the number of species in the AtomVector.
     */
    [[nodiscard]] std::size_t num_species() const noexcept { return m_species_map.size(); }

    /**
     * @brief Convert a species value to the currently used integer to represent that species.
     *
     * @return std::optional<int_fast16_t> Empty if that species is not in the AtomVector.
     */
    [[nodiscard]] std::optional<int_fast16_t> species2z(Symbol const& s) const noexcept {
      for (std::size_t i = 0; i < m_species_map.size(); i++) {
        if (s == m_species_map[i]) {
          return i;
        }
      }
      return {};
    }

    /**
     * @brief Convert the integer used to represent a species to the Symbol used to represent that
     * species.
     *
     * Undefined behaviour if the integer is not used to represent a species in this
     */
    [[nodiscard]] Symbol z2species(int_fast16_t z) const {
      STACK();
      ASSERT(z >= 0 && (std::size_t)z < num_species(), "Species number not in use");
      return m_species_map[z];
    }

    /**
     * @brief Fetch an uninitialised Mat3N large enough to hold the coordinates of every atom
     */
    [[nodiscard]] Mat3N<flt_t> empty_like_z() const {
      return {spatial_dims, static_cast<Eigen::Index>(size())};
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
    [[nodiscard]] Eigen::Map<VecN<int_fast16_t>, Align> z() {
      return {m_z.data(), static_cast<Eigen::Index>(m_z.size())};
    }

    /**
     * @brief Fetch a const Eigen view into the species numbers.
     */
    [[nodiscard]] Eigen::Map<VecN<int_fast16_t> const, Align> z() const {
      return {m_z.data(), static_cast<Eigen::Index>(m_z.size())};
    }

  private:
    std::vector<Symbol> m_species_map;

    std::vector<flt_t, detail::aligned<flt_t, Align>> m_x;
    std::vector<int_fast16_t, detail::aligned<int_fast16_t, Align>> m_z;

  protected:
    friend class bitsery::Access;

    /**
     * @brief Bitsery serialisation
     */
    template <typename S> void serialize(S& s) {
      STACK();
      //
      ASSERT(num_species() < 256, "Too many species");

      s.container(m_species_map, 256);

      ASSERT(size() < 1'000'000, "Too many atoms");

      s.container(m_x, 1'000'000 * 3);
      s.container(m_z, 1'000'000);
    }
  };

}  // namespace otf