#pragma once

#include <fmt/os.h>

#include <optional>

#include "libatom/sim_cell.hpp"

namespace otf::io {

  /**
   * @brief Write the SimCell to a file in extended xyz format and flush the buffer
   *
   * Example:
   *
   * @code{.cpp}
   *
   * #include <fmt/core.h>
   * #include <fmt/os.h>
   *
   * #include "libatom/io/xyz.hpp"
   *
   * auto f = fmt::output_file("dump.xyz");
   *
   * dump_xyz(f, atoms, "A comment!");
   *
   * @endcode
   *
   * @param file File handle
   * @param cell SimCell to take atom data from
   * @param comment Additional comments, must not contain any newline charachters
   */
  void dump_xyz(fmt::ostream& file, SimCell const& cell, std::string_view comment);

}  // namespace otf::io