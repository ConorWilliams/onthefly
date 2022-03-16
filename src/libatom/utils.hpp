#pragma once

#include <string_view>
#include <type_traits>

#include "Eigen/Core"

namespace otf {

  // Compile time constants

#ifndef LIBATOM_SPATIAL_DIMS
  /**
   * @brief Number of spatial dimensions must be 2 or 3
   */
  inline constexpr int spatial_dims = 3;
#else
  /**
   * @brief Number of spatial dimensions must be 2 or 3
   */
  inline constexpr int spatial_dims = LIBATOM_SPATIAL_DIMS;
#  undef LIBATOM_SPATIAL_DIMS
#endif

  // Types aliases

  template <typename T> using Vec = Eigen::Vector<T, spatial_dims>;
  template <typename T> using Mat = Eigen::Matrix<T, spatial_dims, spatial_dims>;

  template <typename T> using VecN = Eigen::Matrix<T, spatial_dims, Eigen::Dynamic>;

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

}  // namespace otf