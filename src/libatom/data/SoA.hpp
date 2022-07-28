#pragma once

#include <type_traits>
#include <utility>

#include <Eigen/Core>

#include "libatom/asserts.hpp"
#include "libatom/detail/atom.hpp"
#include "libatom/utils.hpp"

namespace otf::data {

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
struct AdaptSoA {
public:
  typename Mem::array_t data; ///< Owns its own array.

  /**
   * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
   *
   * For use in ViewSoA hence, models std::vector const-semantics.
   */
  constexpr typename Mem::matrix_ref_t operator()(Mem, int i) {

    ASSERT(i >= 0 && i < data.size() / Mem::size(), "Index out of bounds");

    if constexpr (Mem::is_1x1) {
      return (*data)[i];
    } else {
      return typename Mem::matrix_ref_t{data.data() + i * Mem::size()};
    }
  }

  /**
   * @brief Fetch a view of the ith member stored in the array, tagged dispatch on Mem.
   *
   * For use in ViewSoA hence, models std::vector const-semantics.
   */
  constexpr typename Mem::matrix_cref_t operator()(Mem, int i) const {

    ASSERT(i >= 0 && i < data.size() / Mem::size(), "Index out of bounds");

    if constexpr (Mem::is_1x1) {
      return (*data)[i];
    } else {
      return typename Mem::matrix_cref_t{data.data() + i * Mem::size()};
    }
  }

  /**
   * @brief Fetch a view of the array, tagged dispatch on Mem.
   *
   * For use in ViewSoA hence, models std::vector const-semantics.
   */
  constexpr typename Mem::array_ref_t operator[](Mem) noexcept {
    return data;
  }

  /**
   * @brief Fetch a view of the array, tagged dispatch on Mem.
   *
   * For use in ViewSoA hence, models std::vector const-semantics.
   */
  constexpr typename Mem::array_cref_t operator[](Mem) const noexcept {
    return data;
  }

private:
};

template <typename... Ms>
class SoA : private AdaptSoA<Ms>... {
public:
  /**
   * @brief Construct a new empty SoA
   */
  SoA() = default;

  SoA(SoA&&) = default;

  SoA(SoA const&) = default;

  explicit SoA(int size) : AdaptSoA<Ms>{typename Ms::array_t(size * Ms::size())}..., m_size(size) {}

  template <typename... Mx, std::enable_if_t<!std::is_same_v<SoA, SoA<Mx...>>, bool> = 0>
  explicit SoA(SoA<Mx...> const& other) : AdaptSoA<Ms>{static_cast<AdaptSoA<Ms> const&>(other).data}..., m_size(other.size()) {}

  //   template <typename... Mx, std::enable_if_t<!std::is_same_v<SoA, SoA<Mx...>>, bool> = 0>
  //   explicit SoA(SoA<Mx...>&& other) : AdaptSoA<Ms>(std::move(other.get(Ms{})))..., m_size(other.size()) {}

  using AdaptSoA<Ms>::operator()...;
  using AdaptSoA<Ms>::operator[]...;

  int size() const noexcept {
    return m_size;
  }

  void destructive_resize(int new_size) {
    if (std::exchange(m_size, new_size) != new_size) {
      (static_cast<void>(get(Ms{}).resize(new_size * Ms::size(), Eigen::NoChange)), ...);
    }
  }

private:
  int m_size = 0;

  /* clang-format off */ template <typename...>  friend class SoA; /* clang-format on */
};

} // namespace otf::data