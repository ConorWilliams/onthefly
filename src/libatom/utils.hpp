#pragma once

#include <fmt/chrono.h>
#include <fmt/core.h>

#include <Eigen/Core>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
#include <ratio>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace otf {

// Compile time constants

/**
 * @brief Number of spatial dimensions
 */
inline constexpr std::size_t spatial_dims = 3;

#ifndef LIBATOM_FLOAT_TYPE
#  define LIBATOM_FLOAT_TYPE double
#endif

/**
 * @brief Floating point type used for position, velocity, etc
 */
using floating = LIBATOM_FLOAT_TYPE;

static_assert(std::is_floating_point_v<floating>);

#undef LIBATOM_FLOAT_TYPE

/**
 * @brief The maximum atomic number that any atom can have.
 */
inline constexpr std::size_t max_atomic_num = 111;

// Types aliases

template <typename T>
using Vec3 = Eigen::Array<T, spatial_dims, 1>;
template <typename T>
using Mat3 = Eigen::Matrix<T, spatial_dims, spatial_dims>;

/**
 * @brief Generalised dot-product between to eigen arrays
 */
template <typename E1, typename E2>
auto gdot(Eigen::ArrayBase<E1> const& a, Eigen::ArrayBase<E2> const& b) {
  return (a * b).sum();
}

/**
 * @brief Generic L2-norm squared between two eigen arrays.
 */
template <typename E>
auto norm_sq(Eigen::ArrayBase<E> const& r) {
  return (r * r).sum();
}

/**
 * @brief Generic L2-norm between two eigen arrays.
 */
template <typename E>
auto norm(Eigen::ArrayBase<E> const& r) {
  return std::sqrt(norm_sq(r));
}

/**
 * @brief Remove the soft mode (i.e translating all atoms equally) from the array x.
 */
template <typename E>
auto remove_soft_mode(Eigen::ArrayBase<E> const& x) {
  return x - Eigen::ArrayBase<E>::Ones(x.size()) * x.sum() / x.size();
}

/**
 * @brief Utility to extract first type in a parameter pack.
 */
template <typename T, typename...>
struct First {
  using type = T;
};

/**
 * @brief Compute integer powers of arithmetic types at compile time
 */
template <std::size_t Exp, typename T>
constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> ipow(T base) {
  if constexpr (Exp == 0) {
    return T(1);
  }
  if constexpr (Exp == 1) {
    return base;
  }
  if constexpr (Exp % 2 == 0) {
    return ipow<Exp / 2>(base) * ipow<Exp / 2>(base);
  } else {
    return ipow<Exp - 1>(base) * base;
  }
}

/**
 * @brief A basic implementation of a defer call
 *
 * @code{.cpp}
 *
 * #include <fmt/core.h>
 *
 * #include "libatom/utils.hpp"
 *
 * { // This is a new scope
 *
 *    finally _ = [&] {
 *        fmt::print("Bang!\n");
 *    };
 *
 *    // Do other things here
 *
 * } // <- At scope exit print: Bang!
 *
 * @endcode
 *
 * @tparam F
 */
template <class F>
class [[nodiscard]] finally {
public:
  finally(F f) : m_f(std::move(f)) {}

  finally(const finally&) = delete;
  finally(finally&& other) = delete;
  finally& operator=(const finally&) = delete;
  finally& operator=(finally&&) = delete;

  ~finally() noexcept {
    std::invoke(std::move(m_f));
  }

private:
  F m_f;
};

/**
 * @brief Transparent function wrapper that measures the execution time of a function.
 *
 * The execution time is printed to stdout.
 */
template <typename F, typename... Args>
std::invoke_result_t<F&&, Args&&...> time_call(std::string_view name, F&& f, Args&&... args) {
  //
  auto start = std::chrono::steady_clock::now();

  finally _ = [&] {
    //
    using namespace std::chrono;

    auto elapsed = std::chrono::steady_clock::now() - start;

    auto sec = duration_cast<seconds>(elapsed);

    elapsed -= sec;

    auto mil = duration_cast<milliseconds>(elapsed);

    elapsed -= mil;

    auto mic = duration_cast<microseconds>(elapsed);

    elapsed -= mic;

    auto nan = duration_cast<nanoseconds>(elapsed);

    fmt::print("Timing \"{}\" {:>4} {:>5} {:>5} {:>5}\n", name, sec, mil, mic, nan);
  };

  if constexpr (std::is_void_v<std::invoke_result_t<F&&, Args&&...>>) {
    std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
  } else {
    return std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
  }
}

/**
 * @brief Utility for defining floating chrono types.
 */
template <typename T>
using ftime_t = std::chrono::duration<floating, T>;

/**
 * @brief Quick and dirty timing utility, prints to stdout.
 *
 * @param name Give a name to what you are timing
 * @param f The function you want to time
 * @param timeout The maximum ammount of time you wish to time for, defaults to 1 second
 *
 * @return auto A struct with members .mean and .std sutible for structured binding decomposition
 */
template <typename F>
auto timeit(std::string_view name, F const& f, std::chrono::nanoseconds timeout = std::chrono::seconds{1}) {
  //
  using Duration = std::chrono::high_resolution_clock::duration;

  constexpr std::size_t max_runs = 10'000;

  std::vector<Duration> dt;

  dt.reserve(max_runs);

  Duration elapsed{0};

  fmt::print("Timing \"{}\"...\n", name);

  do {
    auto start = std::chrono::high_resolution_clock::now();

    static_cast<void>(std::invoke(f)); // Discard result

    auto stop = std::chrono::high_resolution_clock::now();

    dt.push_back(stop - start);

    elapsed += stop - start;

  } while (elapsed < timeout && dt.size() < max_runs);

  fmt::print("Performed {} runs in {}\n", dt.size(), elapsed);

  using default_dur = ftime_t<Duration::period>;

  default_dur mean = std::accumulate(dt.begin(), dt.end(), default_dur{0}) / dt.size();

  default_dur nvar{std::accumulate(dt.begin(), dt.end(), 0.0, [&, mean](double acc, default_dur const& val) {
    return acc + (val.count() - mean.count()) * (val.count() - mean.count());
  })};

  default_dur std{std::sqrt(nvar.count() / dt.size())};

  constexpr auto str = "Average time = {} ± {}\n";

  if (mean < ftime_t<std::nano>{1000}) {
    fmt::print(str, ftime_t<std::nano>{mean}, ftime_t<std::nano>{std});
  } else if (mean < ftime_t<std::micro>{1000}) {
    fmt::print(str, ftime_t<std::micro>{mean}, ftime_t<std::micro>{std});
  } else if (mean < ftime_t<std::milli>{1000}) {
    fmt::print(str, ftime_t<std::milli>{mean}, ftime_t<std::milli>{std});
  } else {
    fmt::print(str, ftime_t<std::ratio<1>>{mean}, ftime_t<std::ratio<1>>{std});
  }

  struct Return {
    default_dur mean;
    default_dur std;
  };

  return Return{mean, std};
}

template <typename...>
inline constexpr bool always_false = false;

} // namespace otf
