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
     * @brief Manages a threadsafe stack trace, should not be explicitly instaniciated, use through
     * the corresponding macros. Each StackTrace object is part of an intruisivly linked list which
     * traces the programs call stack.
     */
    class StackTrace {
    public:
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
       * @brief Print the current stacktrace to stderr. In Release this function is a noop
       */
      static void print() {
#ifndef NDEBUG
        _print();
#endif
      }

      /**
       * @brief Construct a new Stack Trace object and append it to the stacktrace's stack
       *
       * @param file name
       * @param func name
       * @param line number
       */
      StackTrace(std::string_view file, std::string_view func, long line) noexcept
          : m_prev{std::exchange(tail()->m_prev, this)}, m_file(file), m_func(func), m_line{line} {}

      /**
       * @brief Destroy the Stack Trace object and pops it from the stack     *
       */
      ~StackTrace() noexcept { tail()->m_prev = m_prev; }

      /**
       * @brief Test two stack traces for equivilance
       */
      bool operator==(StackTrace const& other) {
        return m_line == other.m_line && m_func == other.m_func && m_file == other.m_file;
      }

    private:
      StackTrace* m_prev = nullptr;

      std::string_view m_file{};
      std::string_view m_func{};

      long m_line = 0;

      /**
       * @brief Construct a new Stack Trace object, only for singlton creation
       */
      StackTrace() noexcept = default;

      static void _print();
    };

    /**
     * @brief Kill program printing diagnostics + stacktrace
     *
     * @param expr
     * @param msg
     * @param file
     * @param line
     * @param func
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

#  define OTF_CONCAT_INNER(a, b) a##b
#  define OTF_CONCAT(a, b) OTF_CONCAT_INNER(a, b)

/**
 * @brief RAII like marker, use to store current position on the stack and retrive later through a
 * StackTrace::print()
 */
#  define STACK()                                                  \
    otf::detail::StackTrace OTF_CONCAT(stack_frame, __COUNTER__) { \
      __FILE__, OTF_CURRENT_FUNCTION, __LINE__                     \
    }
/**
 * @brief Use like std c assert but with error message and stacktracing
 */
#  define ASSERT(expr, msg) VERIFY(expr, msg)

#  ifdef eigen_assert
#    error Must not include Eigen before including this header
#  else
#    undef eigen_assert
#    define eigen_assert(x) ASSERT(x, "Eigen embedded assert")
#  endif

#else

#  define STACK() \
    do {          \
    } while (false)

#  define ASSERT(expr, msg) \
    do {                    \
    } while (false)

#endif  // !NDEBUG

}  // namespace otf
