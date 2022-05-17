
#include "libatom/asserts.hpp"

#include <exception>
#include <stdexcept>
#include <string_view>

#include "fmt/core.h"
#include "libatom/utils.hpp"

namespace otf::detail {

  [[noreturn]] void assert_handler(std::string_view expr, std::string_view msg,
                                   std::string_view file, long line, std::string_view func) {
    //
    fmt::print(stderr, "In {}:{}\n", file, line);
    fmt::print(stderr, "From function: {}\n", func);
    fmt::print(stderr, "Assertion \"{}\" failied!\n", expr);
    fmt::print(stderr, "With message: {}\n", msg);

    struct libatom_unrecoverable : std::runtime_error {
      using std::runtime_error::runtime_error;
    };

    throw libatom_unrecoverable("libatom encountered an unrecoverable error");
  }

}  // namespace otf::detail