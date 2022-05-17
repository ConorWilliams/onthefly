#pragma once

#include <cstddef>
#include <cstdint>

#include "libatom/asserts.hpp"
#include "libatom/system/atom_vector.hpp"
#include "libatom/system/member.hpp"
#include "libatom/system/ortho_sim_box.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief A SimCell is a collection of atoms augmented with an OrthoSimBox
   */
  struct SimCell : AtomVector<Position, Frozen, AtomicNum> {
    OrthoSimBox box;
  };

}  // namespace otf