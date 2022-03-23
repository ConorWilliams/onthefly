
#include "libatom/system/simbox.hpp"

#include "libatom/asserts.hpp"
#include "libatom/utils.hpp"

namespace otf {

  OrthoSimBox::OrthoSimBox(Vec<flt_t> const &extents, Vec<bool> const &periodic)
      : m_extents{extents}, m_periodic{periodic}, m_inv_extents(1.0 / extents) {
    //
    STACK();
    VERIFY((m_extents > 0).all(), "OrthoSimBox extents are negative");
  }

}  // namespace otf