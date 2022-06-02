#pragma once

#include <utility>
#include <vector>

#include "libatom/atom.hpp"
#include "libatom/utils.hpp"

namespace otf::env {

  /**
   * @brief A geometry models a local distribution of Atoms.
   *
   * The distribution is centred on the first atom.
   */
  template <typename... Mems> class Geometry : public AtomVector<Mems...> {};

}  // namespace otf::env