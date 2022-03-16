
#include "asserts.hpp"

#include <string_view>

#include "fmt/core.h"

namespace otf::detail {

  void StackTrace::_print() {
    //

    fmt::print(stderr, "Stacktrace:\n");

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

    int c = 1;
    int r = 0;

    for (StackTrace* f = tail()->prev; f != nullptr; f = f->prev) {
      //
      if (f->prev && *f == *(f->prev)) {
        r += 1;
        c += 1;
      } else if (r > 0) {
        fmt::print(stderr, "  #{:<4} {:{}} {:5} {}\n", c++ - r, f->_file, len, f->_line, f->_func);
        fmt::print(stderr, "|                                 |\n");
        fmt::print(stderr, "| {:>4} layer(s) of recursion      |\n", r);
        fmt::print(stderr, "|                                 |\n");
        fmt::print(stderr, "  #{:<4} {:{}} {:5} {}\n", c++, f->_file, len, f->_line, f->_func);
        r = 0;
      } else {
        fmt::print(stderr, "  #{:<4} {:{}} {:5} {}\n", c++, f->_file, len, f->_line, f->_func);
      }
    }
  }

  [[noreturn]] void assert_handler(std::string_view expr, std::string_view msg,
                                   std::string_view file, long line, std::string_view func) {
    //
    fmt::print(stderr, "In {}:{} from function: {}\n", file, line, func);
    fmt::print(stderr, "Assertion \"{}\" failied with message: {}\n", expr, msg);

    otf::detail::StackTrace::print();

    std::terminate();
  }

}  // namespace otf::detail