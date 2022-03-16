
#include "asserts.hpp"

#include <string_view>

#include "fmt/core.h"

namespace otf {

  void StackTrace::print() {
    //
    int count = 1;

    fmt::print(stderr, "\nStacktrace:\n");

    // Find longest common prefix

    std::string_view prefix = tail()->prev ? tail()->prev->_file : "";

    for (StackTrace* f = tail()->prev; f != nullptr; f = f->prev) {
      prefix = common_prefix(prefix, f->_file);
    }

    // Trim all filenames and compute length of the longest

    std::size_t len = 0;

    for (StackTrace* f = tail()->prev; f != nullptr; f = f->prev) {
      f->_file.remove_prefix(prefix.size());

      len = std::max(len, f->_file.size());
    }

    // Print stacktrace

    for (StackTrace* f = tail()->prev; f != nullptr; f = f->prev) {
      fmt::print(stderr, "  #{:<4} {:{}} {:5} {}\n", count++, f->_file, len, f->_line, f->_func);
    }
  }

}  // namespace otf