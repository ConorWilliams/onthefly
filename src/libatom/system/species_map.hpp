#pragma once

#include <algorithm>
#include <array>
#include <cstddef>  // size_t
#include <cstdint>  // fast16_t
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "libatom/asserts.hpp"

namespace otf {

  /**
   * @brief Class to manage the mapping of string-symbols <-> integer-keys.
   *
   * The atom integer-keys's ("z"s) are constructed sequentially such that each new atom type gets
   * a new integer-keys. The species symbol can be reconstructed using .z2species().
   */
  class SpeciesMap {
  public:
    /**
     * @brief Get the number of species in the SpeciesMap.
     */
    [[nodiscard]] std::size_t size() const noexcept { return m_species_map.size(); }

    /**
     * @brief Convert a Species to the currently used integer to represent that species.
     *
     * @return std::optional<std::size_t> Empty if that species is not in the SpeciesMap.
     */
    [[nodiscard]] std::optional<std::size_t> species2z(std::string_view s) const noexcept {
      for (std::size_t i = 0; i < m_species_map.size(); i++) {
        if (s == m_species_map[i]) {
          return i;
        }
      }
      return {};
    }

    /**
     * @brief Same as .species2z() but insert species if not currently in the map.
     */
    std::size_t species2z_or_insert(std::string_view const& s) {
      if (std::optional z = species2z(s)) {
        return *z;
      } else {
        m_species_map.emplace_back(s);
        return size() - 1;
      }
    }

    /**
     * @brief Convert the integer used to represent a species to the Symbol used to represent that
     * species.
     *
     * Undefined behaviour if the integer is not used to represent a species in this SpeciesMap.
     */
    [[nodiscard]] std::string_view z2species(std::size_t z) const {
      ASSERT(z < size(), "Species number not in use");
      return m_species_map[z];
    }

  private:
    std::vector<std::string> m_species_map;
  };

}  // namespace otf
