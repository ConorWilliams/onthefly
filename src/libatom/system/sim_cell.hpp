#pragma once

#include <cstddef>
#include <cstdint>

#include "libatom/asserts.hpp"
#include "libatom/system/atom_array.hpp"
#include "libatom/system/member.hpp"
#include "libatom/system/ortho_sim_box.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief A SimCell is a collection of atoms augmented with an OrthoSimBox
   */
  class SimCell : public AtomArray<Position, Frozen, AtomicNum, Gradient> {
  public:
    OrthoSimBox box;

    SimCell(OrthoSimBox const& arg) : box{arg} {}

  private:
  };

}  // namespace otf