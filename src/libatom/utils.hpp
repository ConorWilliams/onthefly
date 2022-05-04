#pragma once

#include <chrono>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "fmt/chrono.h"
#include "fmt/core.h"

//

#include "libatom/asserts.hpp"
// ^ must procede Eigen
#include "Eigen/Core"

namespace otf {

  // Compile time constants

#ifndef LIBATOM_SPATIAL_DIMS
#  define LIBATOM_SPATIAL_DIMS 3
#endif

  /**
   * @brief Number of spatial dimensions must be 2 or 3
   */
  inline constexpr int spatial_dims{LIBATOM_SPATIAL_DIMS};

  static_assert(spatial_dims == 2 || spatial_dims == 3, "Invalid number of dimensions");

#undef LIBATOM_SPATIAL_DIMS

#ifndef LIBATOM_FLOAT_TYPE
#  define LIBATOM_FLOAT_TYPE double
#endif

  /**
   * @brief Floating point type used for position, velocity, etc
   */
  using flt_t = LIBATOM_FLOAT_TYPE;

  static_assert(std::is_floating_point_v<flt_t>);

#undef LIBATOM_FLOAT_TYPE

  // Types aliases

  template <typename T> using Vec = Eigen::Array<T, spatial_dims, 1>;
  template <typename T> using Mat = Eigen::Array<T, spatial_dims, spatial_dims>;

  template <typename T> using VecN = Eigen::Array<T, Eigen::Dynamic, 1>;
  template <typename T> using Mat3N = Eigen::Array<T, spatial_dims, Eigen::Dynamic>;

  /**
   * @brief Return the longest prefix of the two inputs strings.
   *
   * @param Input string a
   * @param Input string b
   *
   * @return A substring of the input string_view "a" containing the longest common prefix of a
   * and b.
   */
  std::string_view common_prefix(std::string_view a, std::string_view b);

  /**
   * @brief Generic L2-norm squared between two eigen arrays.
   */
  template <typename E> double norm_sq(Eigen::ArrayBase<E> const& r) { return (r * r).sum(); }

  /**
   * @brief Generic L2-norm between two eigen arrays.
   */
  template <typename E> double norm(Eigen::ArrayBase<E> const& r) { return std::sqrt(norm_sq(r)); }

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
   * @brief Quick and dirty timing utility, prints to stdout.
   *
   * @param name Give a name to what you are timing
   * @param f The function you want to time
   * @param timeout The maximum ammount of time you wish to time for, defaults to 1 second
   *
   * @return auto A struct with members .mean and .std sutible for structured binding decomposition
   */
  template <typename F> auto timeit(std::string_view name, F const& f,
                                    std::chrono::seconds timeout = std::chrono::seconds{1}) {
    //
    using Duration = std::chrono::high_resolution_clock::duration;

    constexpr std::size_t max_runs = 10'000;

    std::vector<Duration> dt;

    dt.reserve(max_runs);

    Duration elapsed{0};

    fmt::print("Timing \"{}\"...\n", name);

    do {
      auto start = std::chrono::high_resolution_clock::now();

      static_cast<void>(std::invoke(f));  // Discard result

      auto stop = std::chrono::high_resolution_clock::now();

      dt.push_back(stop - start);

      elapsed += stop - start;

    } while (elapsed < timeout && dt.size() < max_runs);

    fmt::print("Performed {} runs in {}\n", dt.size(), elapsed);

    using floating_ns = std::chrono::duration<double, Duration::period>;

    floating_ns mean = std::accumulate(dt.begin(), dt.end(), floating_ns{0}) / dt.size();

    floating_ns nvar{
        std::accumulate(dt.begin(), dt.end(), 0.0, [&, mean](double acc, floating_ns const& val) {
          return acc + (val.count() - mean.count()) * (val.count() - mean.count());
        })};

    floating_ns std{std::sqrt(nvar.count() / dt.size())};

    fmt::print("Average time = {} ± {}\n", mean, std);

    struct Return {
      floating_ns mean;
      floating_ns std;
    };

    return Return{mean, std};
  }

}  // namespace otf

#include "bitsery/bitsery.h"
#include "bitsery/brief_syntax.h"

namespace bitsery {
  /**
   * @brief Serializeation of Vec<T> using bitsery
   */
  template <typename S, typename T> void serialize(S& s, otf::Vec<T>& v) {
    if constexpr (otf::spatial_dims == 3) {
      s(v[0], v[1], v[2]);
    } else if constexpr (otf::spatial_dims == 2) {
      s(v[0], v[1]);
    } else {
      // Fallback to loop
      for (auto&& elem : v) {
        s(elem);
      }
    }
  }
}  // namespace bitsery