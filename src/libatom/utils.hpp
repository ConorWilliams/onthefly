#pragma once

#include <cmath>
#include <string_view>
#include <type_traits>

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

#undef LIBATOM_FLOAT_TYPE

#ifndef LIBATOM_FLOAT_TYPE
#  define LIBATOM_FLOAT_TYPE double
#endif

  /**
   * @brief Floating point type used for position, velocity, etc
   */
  using float_t = LIBATOM_FLOAT_TYPE;

  static_assert(std::is_floating_point_v<float_t>);

#ifndef LIBATOM_FLOAT_ACC_TYPE
#  define LIBATOM_FLOAT_ACC_TYPE LIBATOM_FLOAT_TYPE
#endif

  /**
   * @brief Floating point type used for accumulating float_t
   */
  using float_acc_t = LIBATOM_FLOAT_ACC_TYPE;

  static_assert(std::is_floating_point_v<float_acc_t>);

#undef LIBATOM_FLOAT_TYPE
#undef LIBATOM_FLOAT_ACC_TYPE

  // Types aliases

  template <typename T> using Vec = Eigen::Array<T, spatial_dims, 1>;
  template <typename T> using Mat = Eigen::Array<T, spatial_dims, spatial_dims>;

  template <typename T> using VecN = Eigen::Array<T, spatial_dims, Eigen::Dynamic>;

  /**
   * @brief A C++23 version of std::exchange with constexpr+noexcept
   *
   * @param obj Reference to object that will be replaced.
   * @param new_value To write to obj.
   * @return Old value of obj.
   */
  template <typename T, typename U = T> constexpr T exchange(T& obj, U&& new_value) noexcept(
      std::is_nothrow_move_constructible_v<T>&& std::is_nothrow_assignable_v<T&, U>) {
    T old_value = std::move(obj);
    obj = std::forward<U>(new_value);
    return old_value;
  }

  /**
   * @brief Return the longest prefix of the two inputs strings
   *
   * @param Input string a
   * @param Input string b
   *
   * @return A view substring of the input string_view a containing the longest common prefix of a
   * and b.
   */
  std::string_view common_prefix(std::string_view a, std::string_view b);

  /**
   * @brief Generic L2-norm squared between two eigen arrays
   */
  template <typename E> double norm_sq(Eigen::ArrayBase<E> const& r) { return (r * r).sum(); }

  /**
   * @brief Generic L2-norm between two eigen arrays
   */
  template <typename E> double norm(Eigen::ArrayBase<E> const& r) { return std::sqrt(norm_sq(r)); }

}  // namespace otf