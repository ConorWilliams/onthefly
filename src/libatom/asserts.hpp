#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>

#include "fmt/core.h"
#include "libatom/external/current_function.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief Manages a threadsafe stack trace, should not be explicitly instaniciated, use through
   * the corresponding macros. Each StackTrace object is part of an intruisivly linked list which
   * traces the stack.
   */
  struct StackTrace {
    //
    StackTrace* prev = nullptr;

    std::string_view _file{};
    std::string_view _func{};

    long _line = 0;

    /**
     * @brief Per-thread singleton which acts as the tail of the list
     *
     * @return StackTrace* singleton instance
     */
    static StackTrace* tail() {
      thread_local StackTrace tail{};
      return &tail;
    }

    /**
     * @brief Print the current stacktrace to stderr.
     */
    static void print();

    /**
     * @brief Construct a new Stack Trace object and append it to the stacktrace's stack
     *
     * @param file name
     * @param func name
     * @param line number
     */
    StackTrace(std::string_view file, std::string_view func, long line) noexcept
        : prev{otf::exchange(tail()->prev, this)}, _file(file), _func(func), _line{line} {}

    /**
     * @brief Destroy the Stack Trace object and pops it from the stack     *
     */
    ~StackTrace() noexcept { tail()->prev = prev; }

  private:
    /**
     * @brief Construct a new Stack Trace object, only for singlton creation
     */
    StackTrace() noexcept = default;
  };

#ifndef NDEBUG

#  define OTF_CONCAT_INNER(a, b) a##b
#  define OTF_CONCAT(a, b) OTF_CONCAT_INNER(a, b)

/**
 * @brief RAII like marker, use to store current position on the stack and retrive later through a
 * stack trace *
 */
#  define STACK()                                          \
    otf::StackTrace OTF_CONCAT(stack_frame, __COUNTER__) { \
      __FILE__, BOOST_CURRENT_FUNCTION, __LINE__           \
    }

#  define ASSERT(expr, msg)                                                               \
    do {                                                                                  \
      if (!(expr)) {                                                                      \
        fmt::print(stderr, "\nAssertion \"{}\" failied with message: {}\n", #expr, #msg); \
        fmt::print(stderr, "{}: on line {} in function: {}\n", __FILE__, __LINE__,        \
                   BOOST_CURRENT_FUNCTION);                                               \
        otf::StackTrace::print();                                                         \
        std::terminate();                                                                 \
      }                                                                                   \
    } while (false)

#else

#  define OTF_STACK() \
    do {              \
    } while (false)

#  define OTF_CHECK(expr) \
    do {                  \
    } while (false)

#endif  // !NDEBUG

}  // namespace otf
