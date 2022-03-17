
#include "libatom/asserts.hpp"

#include <exception>
#include <string_view>

#include "fmt/core.h"

namespace otf::detail {

  void StackTrace::_print() {
    //

    fmt::print(stderr, "Stacktrace:\n");

    // Find longest common prefix

    std::string_view prefix = tail()->m_prev ? tail()->m_prev->m_file : "";

    for (StackTrace* f = tail()->m_prev; f != nullptr; f = f->m_prev) {
      prefix = common_prefix(prefix, f->m_file);
    }

    // Trim all filenames and compute length of the longest

    std::size_t len = 0;

    for (StackTrace* f = tail()->m_prev; f != nullptr; f = f->m_prev) {
      f->m_file.remove_prefix(prefix.size());

      len = std::max(len, f->m_file.size());
    }

    // Print stacktrace

    int c = 1;
    int r = 0;

    for (StackTrace* f = tail()->m_prev; f != nullptr; f = f->m_prev) {
      //
      if (f->m_prev && *f == *(f->m_prev)) {
        r += 1;
        c += 1;
      } else if (r > 0) {
        fmt::print(stderr, "  #{:<4} {:{}} {:5} {}\n", c++ - r, f->m_file, len, f->m_line,
                   f->m_func);
        fmt::print(stderr, "|                                 |\n");
        fmt::print(stderr, "| {:>4} layer(s) of recursion      |\n", r);
        fmt::print(stderr, "|                                 |\n");
        fmt::print(stderr, "  #{:<4} {:{}} {:5} {}\n", c++, f->m_file, len, f->m_line, f->m_func);
        r = 0;
      } else {
        fmt::print(stderr, "  #{:<4} {:{}} {:5} {}\n", c++, f->m_file, len, f->m_line, f->m_func);
      }
    }
  }

  struct libatom_unrecoverable : std::exception {};

  [[noreturn]] void assert_handler(std::string_view expr, std::string_view msg,
                                   std::string_view file, long line, std::string_view func) {
    //
    fmt::print(stderr, "In {}:{} from function: {}\n", file, line, func);
    fmt::print(stderr, "Assertion \"{}\" failied with message: {}\n", expr, msg);

    otf::detail::StackTrace::print();

    throw libatom_unrecoverable{};
  }

}  // namespace otf::detail