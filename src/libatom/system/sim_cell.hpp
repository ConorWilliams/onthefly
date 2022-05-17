#pragma once

#include <cstddef>
#include <cstdint>

#include "libatom/asserts.hpp"
#include "libatom/system/atom_vector.hpp"
#include "libatom/system/ortho_sim_box.hpp"
#include "libatom/system/species_map.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief A SimCell is a collection of  atoms, a species map and an OrthoSimBox
   */
  class SimCell {
  public:
    SpeciesMap map;

    OrthoSimBox box;

    AtomVector active;
    AtomVector frozen;

    /**
     * @brief Construct a new Sim Cell object containing no atoms
     */
    explicit SimCell(OrthoSimBox const& box) : box{box} {}

    /**
     * @brief Get the total number of atoms in active + frozen
     */
    std::size_t size() const noexcept { return active.size() + frozen.size(); }

    // VecN<double> active_disp(VecN<double> const &others) const;
  };

}  // namespace otf