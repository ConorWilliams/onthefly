#pragma once

#include <cstddef>
#include <type_traits>
#include <utility>

#include <Eigen/Core>

#include "libatom/asserts.hpp"
#include "libatom/detail/atom.hpp"
#include "libatom/utils.hpp"

namespace otf::data2 {

/**
 * @brief A base type to derive from for defining members of an Atom type.
 *
 * Members must be matrices of arithmetic types or default constructible 1x1 matricies.
 *
 * 1x1 matricies are unwrapped into scalars
 *
 * @tparam Scalar This member represents a matrix of Scalar elements.
 * @tparam Rows Number of rows in this member.
 * @tparam Cols Number of colums in this member.
 * @tparam Rep The Eigen3 template, Eigen::[matrix||array], to use for this member.
 */
template <typename Scalar, int Rows = 1, int Cols = 1, template <typename, auto...> typename Rep = Eigen::Matrix>
struct MemTag {
  /** @brief True if this member represents a 1x1 matrix. */
  static constexpr bool is_1x1 = Rows == 1 && Cols == 1;

  static_assert(std::is_arithmetic_v<Scalar> || (Rows == 1 && Cols == 1), "Non-scalar members must be arithmetic.");

  static_assert(std::is_default_constructible_v<Scalar>, "Scalar members must be default constructable.");

  static_assert(Rows > 0 && Cols > 0, "Invalid member extents.");

  /** @brief This member represents a matrix of elements of scalar_t. */
  using scalar_t = Scalar;

  /** @brief The matrix type that this member represents (1x1 matrices are unwrapped to scalars). */
  using matrix_t = std::conditional_t<is_1x1, Scalar, Rep<Scalar, Rows, Cols>>;
  /** @brief A reference-like type to a matrix_t. */
  using matrix_ref_t = std::conditional_t<is_1x1, Scalar&, Eigen::Map<matrix_t>>;
  /** @brief A const-reference-like type to a const matrix_t. */
  using matrix_cref_t = std::conditional_t<is_1x1, Scalar const&, Eigen::Map<matrix_t const>>;

  /** @brief The Eigen type used to store a dynamic collection of contiguous matrix_t. */
  using array_t = Eigen::Array<Scalar, Eigen::Dynamic, 1>;
  /** @brief A reference-like type to the underlying array_t. */
  using array_ref_t = Eigen::ArrayBase<array_t>&;
  /** @brief A const-reference-like type to the underlying array_t. */
  using array_cref_t = Eigen::ArrayBase<array_t> const&;

  /**
   * @brief Get the number of elements in the matrix_t.
   */
  static constexpr int size() {
    return Rows * Cols;
  }
};

template <typename Mem>
struct Adaptor {
public:
  Adaptor() = default;

  Adaptor(Adaptor&&) = default;

  Adaptor(Adaptor const&) = default;

  /**
   * Possibilities
   *
   *  void foo1(Adaptor<Pos>)
   *  void foo2(Adaptor<Pos> &)
   *  void foo3(Adaptor<Pos> const &)
   *
   *  foo1(Adaptor<Pos>{})          use defaulted copy/move
   *  foo1(Adaptor<Pos &>{})        ok to do implicitly as user wants a copy
   *  foo1(Adaptor<Pos const&>{})   ok to do implicitly as user wants a copy
   *
   *  foo2(Adaptor<Pos>{})          use defaulted copy/move, compiler will not allow ref to temp
   *  foo2(Adaptor<Pos &>{})        compiler will not allow ref to temp
   *  foo2(Adaptor<Pos const&>{})   compiler will not allow ref to temp
   *
   *  foo3(Adaptor<Pos>{})          Just a const ref
   *  foo3(Adaptor<Pos &>{})        must not allow implicit construction from a view else this would create a temporary
   *  foo3(Adaptor<Pos const&>{})   same as above ^
   *
   */

  // If constructing from a const view then deference their pointer. Pass views by value.
  explicit Adaptor(Adaptor<Mem const&> other) : m_data(*other.m_data_ptr) {}

  // If constructing from a view then deference their pointer. Pass views by value.
  explicit Adaptor(Adaptor<Mem&> other) : m_data(*other.m_data_ptr) {}

  // Owning specific.
  explicit Adaptor(std::size_t size) : m_data(size * Mem::size()) {}

  Adaptor& operator=(Adaptor const&) = default;

  Adaptor& operator=(Adaptor&&) = default;

  Adaptor& operator=(Adaptor<Mem const&> other) {
    m_data = *other.m_data;
    return *this;
  }

  Adaptor& operator=(Adaptor<Mem&> other) {
    m_data = *other.m_data;
    return *this;
  }

  /**
   * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
   *
   * This is an owning Adaptor hence, model value const-semantics.
   */
  constexpr typename Mem::matrix_ref_t operator()(Mem, int i) {

    ASSERT(i >= 0 && i < m_data.size() / Mem::size(), "Index out of bounds");

    if constexpr (Mem::is_1x1) {
      return m_data[i];
    } else {
      return typename Mem::matrix_ref_t{m_data.data() + i * Mem::size()};
    }
  }

  /**
   * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
   *
   * This is an owning Adaptor hence, model value const-semantics
   */
  constexpr typename Mem::matrix_cref_t operator()(Mem, int i) const {

    ASSERT(i >= 0 && i < m_data.size() / Mem::size(), "Index out of bounds");

    if constexpr (Mem::is_1x1) {
      return m_data[i];
    } else {
      return typename Mem::matrix_cref_t{m_data.data() + i * Mem::size()};
    }
  }

  /**
   * @brief Fetch a view of the array, tagged dispatch on Mem.
   *
   * This is an owning Adaptor hence, model value const-semantics
   */
  constexpr typename Mem::array_ref_t operator[](Mem) noexcept {
    return m_data;
  }

  /**
   * @brief Fetch a view of the array, tagged dispatch on Mem.
   *
   * This is an owning Adaptor hence, model value const-semantics
   */
  constexpr typename Mem::array_cref_t operator[](Mem) const noexcept {
    return m_data;
  }

private:
  typename Mem::array_t m_data; ///< Owns its own array.

  friend struct Adaptor<Mem&>;
  friend struct Adaptor<Mem const&>;
};

template <typename Mem>
struct Adaptor<Mem&> {
public:
  Adaptor() = default;

  Adaptor(Adaptor&&) = default;

  Adaptor(Adaptor const&) = default;

  // If constructing from am owning Adaptor then take address of their m_data member, no explicit for construction of a view.
  Adaptor(Adaptor<Mem>& other) : m_data_ptr(&other.m_data) {}

  // Cannot construct from a const view.
  Adaptor(Adaptor<Mem const&> other) = delete;

  Adaptor& operator=(Adaptor const&) = default;

  Adaptor& operator=(Adaptor&&) = default;

  // Cannot assign to a const view
  Adaptor& operator=(Adaptor<Mem const&> other) = delete;

  // Can only assign to a reference to an owning view
  Adaptor& operator=(Adaptor<Mem>& other) {
    m_data_ptr = &other.m_data;
    return *this;
  }

  /**
   * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
   *
   * This is not an owning Adaptor hence, model pointer const-semantics.
   */
  constexpr typename Mem::matrix_ref_t operator()(Mem, int i) const {

    ASSERT(i >= 0 && i < m_data_ptr->size() / Mem::size(), "Index out of bounds");

    if constexpr (Mem::is_1x1) {
      return (*m_data_ptr)[i];
    } else {
      return typename Mem::matrix_ref_t{m_data_ptr->m_data() + i * Mem::size()};
    }
  }

  /**
   * @brief Fetch a view of the array, tagged dispatch on Mem.
   *
   * This is not an owning Adaptor hence, model pointer const-semantics.
   */
  constexpr typename Mem::array_ref_t operator[](Mem) const noexcept {
    return *m_data_ptr;
  }

private:
  Eigen::ArrayBase<typename Mem::array_t>* m_data_ptr; ///< Pointer to the viewed array.

  friend struct Adaptor<Mem>;
  friend struct Adaptor<Mem const&>;
};

template <typename Mem>
struct Adaptor<Mem const&> {
public:
  Adaptor() = default;

  Adaptor(Adaptor&&) = default;

  Adaptor(Adaptor const&) = default;

  // If constructing from a non const view then just copy the pointer, no explicit for construction of a view.
  Adaptor(Adaptor<Mem&> other) : m_data_ptr(other.m_data_ptr) {}

  // If constructing from am owning Adaptor then take address of their m_data member, no explicit for construction of a view.
  Adaptor(Adaptor<Mem> const& other) : m_data_ptr(&other.m_data) {}

  Adaptor& operator=(Adaptor const&) = default;

  Adaptor& operator=(Adaptor&&) = default;

  Adaptor& operator=(Adaptor<Mem&> other) {
    m_data_ptr = other.m_data_ptr;
    return *this;
  }
  // Can only assign to a reference to an owning view
  Adaptor& operator=(Adaptor<Mem> const& other) {
    m_data_ptr = &other.m_data;
    return *this;
  }

  /**
   * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
   *
   * This is not an owning Adaptor hence, model const-pointer const-semantics.
   */
  constexpr typename Mem::matrix_cref_t operator()(Mem, int i) const {

    ASSERT(i >= 0 && i < m_data_ptr->size() / Mem::size(), "Index out of bounds");

    if constexpr (Mem::is_1x1) {
      return (*m_data_ptr)[i];
    } else {
      return typename Mem::matrix_cref_t{m_data_ptr->m_data() + i * Mem::size()};
    }
  }

  /**
   * @brief Fetch a view of the array, tagged dispatch on Mem.
   *
   * This is not an owning Adaptor hence, model const-pointer const-semantics.
   */
  constexpr typename Mem::array_cref_t operator[](Mem) const noexcept {
    return *m_data_ptr;
  }

private:
  Eigen::ArrayBase<typename Mem::array_t> const* m_data_ptr; ///< Pointer to the viewed array.

  friend struct Adaptor<Mem>;
  friend struct Adaptor<Mem&>;
};

template <typename Mem>
struct Adaptor<Mem&&> {
  static_assert(always_false<Mem>, "Rvalue reference to member is illegal.");
};

template <typename Mem>
struct Adaptor<Mem const&&> {
  static_assert(always_false<Mem>, "Rvalue reference to member is illegal.");
};

template <typename T>
using remove_cref_t = std::remove_const_t<std::remove_reference_t<T>>;

template <typename... Ms>
class SoA;

namespace detail {

template <typename, typename>
struct different_SoA : std::false_type {};

template <typename... Ms>
struct different_SoA<SoA<Ms...>, SoA<Ms...>> : std::false_type {};

template <typename... Ms, typename... Mx>
struct different_SoA<SoA<Ms...>, SoA<Mx...>> : std::true_type {};

} // namespace detail

template <typename... Ms>
class SoA : private Adaptor<Ms>... {
private:
  static constexpr bool owns_none = (std::is_reference_v<Ms> && ...);
  static constexpr bool owns_all = (!std::is_reference_v<Ms> && ...);

  template <typename T>
  static constexpr bool different_SoA_v = detail::different_SoA<SoA<Ms...>, remove_cref_t<T>>::value;

public:
  /**
   * @brief Construct a new empty SoA
   */
  SoA() = default;

  SoA(SoA&&) = default;

  SoA(SoA const&) = default;

  /**
   * @brief Construct a new SoA conatining 'size' default initializes atoms.
   *
   * Only SFINE enabled if this SoA owns all its arrays.
   */
  template <bool Owning = owns_all>
  explicit SoA(int size, std::enable_if_t<Owning>* = 0) : Adaptor<Ms>(size)..., m_size(size) {}

  /**
   * @brief Implicitly construct a new SoA object from SoA 'other' with different members.
   *
   * Only SFINE enabled if this SoA owns non of its arrays.
   *
   * .. note::
   *    The implementation may ``std::move`` ``other`` multiple times but this is ok as the Adaptor constructor will only move its
   *    corresponding base slice.
   *
   */
  template <typename T, typename = std::enable_if_t<different_SoA_v<T> && owns_none>>
  SoA(T&& other) : Adaptor<Ms>(std::forward<T>(other))..., m_size(other.size()) {}

  /**
   * @brief Explicitly construct a new SoA object from SoA 'other' with different members.
   *
   * SFINE enabled if this SoA owns some of its arrays.
   *
   * .. note::
   *    The implementation may ``std::move`` ``other`` multiple times but this is ok as the Adaptor constructor will only move its
   *    corresponding base slice.
   */
  template <typename T, typename = std::enable_if_t<different_SoA_v<T>>, typename = void>
  explicit SoA(T&& other) : Adaptor<Ms>(std::forward<T>(other))..., m_size(other.size()) {}

  //

  using Adaptor<Ms>::operator()...;
  using Adaptor<Ms>::operator[]...;

  int size() const noexcept {
    return m_size;
  }
  template <bool Owning = owns_all>
  void destructive_resize(int new_size, std::enable_if_t<Owning>* = 0) {
    if (std::exchange(m_size, new_size) != new_size) {
      (static_cast<void>(get(Ms{}).resize(new_size * Ms::size(), Eigen::NoChange)), ...);
    }
  }

private:
  int m_size = 0;

  /* clang-format off */ template <typename...>  friend class SoA; /* clang-format on */
};

} // namespace otf::data2