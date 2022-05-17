#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>

#include "libatom/external/current_function.hpp"

namespace otf {

  namespace detail {
    /**
     * @brief Kill program printing diagnostics + stacktrace, should not be inlined to reduce code
     * size
     */
    [[noreturn]] void assert_handler(std::string_view expr, std::string_view msg,
                                     std::string_view file, long line, std::string_view func);
  }  // namespace detail

// Like ASSERT but on in release build
#define VERIFY(expr, msg)                                                                 \
  do {                                                                                    \
    if (!(expr)) {                                                                        \
      otf::detail::assert_handler(#expr, #msg, __FILE__, __LINE__, OTF_CURRENT_FUNCTION); \
    }                                                                                     \
  } while (false)

#ifndef NDEBUG

/**
 * @brief Use like std c assert but with error message and stacktracing
 */
#  define ASSERT(expr, msg) VERIFY(expr, msg)

#else

#  define ASSERT(expr, msg) \
    do {                    \
    } while (false)

#endif  // !NDEBUG

}  // namespace otf
