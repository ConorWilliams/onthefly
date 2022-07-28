#pragma once

#include <cstddef>
#include <type_traits>
#include <utility>

#include <Eigen/Core>

#include "libatom/asserts.hpp"

#include "libatom/data/SoA.hpp"
#include "libatom/utils.hpp"

namespace otf::data {

/**
 * @brief An Aggregate class that manages a view into a single array.
 *
 * @tparam Mem Tag class (may be const qualified), identifies member type stored in array.
 */
template <typename Mem>
struct AdaptViewSoA {
private:
  using plain_Mem = std::remove_const_t<Mem>;
  // Get Eigen base type used to store Mem.
  using base_t = Eigen::ArrayBase<typename plain_Mem::array_t>;
  // Correctly const-qualified base_t.
  using view_t = std::conditional_t<std::is_const_v<Mem>, base_t const, base_t>;

  // Conditionally const-qualified matrix_ref_t.
  using matrix_ref_t = std::conditional_t<std::is_const_v<Mem>, typename plain_Mem::matrix_cref_t, typename plain_Mem::matrix_ref_t>;
  // Conditionally const-qualified array_ref_t.
  using array_ref_t = std::conditional_t<std::is_const_v<Mem>, typename plain_Mem::array_cref_t, typename plain_Mem::array_ref_t>;

public:
  view_t* data_ptr; ///< Pointer to the viewed array.

  /**
   * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
   *
   * For use in ViewSoA hence, models std::span const-semantics.
   */
  constexpr matrix_ref_t operator()(plain_Mem, int i) const {

    ASSERT(i >= 0 && i < data_ptr->size() / Mem::size(), "Index out of bounds");

    if constexpr (plain_Mem::is_1x1) {
      return (*data_ptr)[i];
    } else {
      return matrix_ref_t{data_ptr->data() + i * plain_Mem::size()};
    }
  }

  /**
   * @brief Fetch a view of the array, tagged dispatch on Mem.
   *
   * For use in ViewSoA hence, models std::span const-semantics.
   */
  constexpr array_ref_t operator[](plain_Mem) const noexcept {
    return *data_ptr;
  }
};

template <typename... Ms>
class ViewSoA : private AdaptViewSoA<Ms>... {
public:
  /**
   * @brief Construct a new empty ViewSoA
   */
  ViewSoA() = default;

  ViewSoA(ViewSoA&&) = default;

  ViewSoA(ViewSoA const&) = default;

  template <typename... Mx>
  ViewSoA(SoA<Mx...>& other) : AdaptViewSoA<Ms>{&other[std::remove_const_t<Ms>{}]}..., m_size(other.size()) {}

  int size() const noexcept {
    return m_size;
  }

  using AdaptViewSoA<Ms>::operator()...;
  using AdaptViewSoA<Ms>::operator[]...;

private:
  int m_size = 0;
};

} // namespace otf::data