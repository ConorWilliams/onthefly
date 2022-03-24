#pragma once

#include <optional>

#include "fmt/os.h"
#include "libatom/system/simcell.hpp"

namespace otf {

  /**
   * @brief Write the SimCell to a file in extended xyz format and flush
   *
   * @param file File handle
   * @param cell SimCell to take atom data from
   * @param time Optional time parameter
   */
  void to_xyz(fmt::ostream& file, SimCell const& cell, std::optional<flt_t> time = std::nullopt);

}  // namespace otf