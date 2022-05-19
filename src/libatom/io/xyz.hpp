#pragma once

#include <optional>

#include "fmt/os.h"
#include "libatom/system/sim_cell.hpp"

namespace otf {

  /**
   * @brief Write the SimCell to a file in extended xyz format and flush the buffer
   *
   * Example:
   *
   * @code{.cpp}
   *
   * #include <fmt/core.h>
   * #include <fmt/ostream.h>
   *
   * #include "libatom/io/xyz.hpp"
   *
   * dump_xyz(fmt::output_file("dump.xyz"), atoms, "A comment!");
   *
   * @endcode
   *
   * @param file File handle
   * @param cell SimCell to take atom data from
   * @param time Additional comments, must not contain any newline charachters
   */
  void dump_xyz(fmt::ostream& file, SimCell const& cell, std::string_view comment);

}  // namespace otf