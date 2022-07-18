#pragma once

#include <cstddef>
#include <type_traits>
#include <utility>

#include "libatom/utils.hpp"

/**
 * @brief Internal namespace, not part of the API.
 */
namespace otf::detail {

  /**
   * @brief Stores a single member of an Atom and allows tagged access.
   *
   * This is an internal type used by Atom.
   *
   * @tparam Tag An empty type, derived from otf::MemTag.
   */
  template <typename Tag> struct AtomMem {
    static_assert(std::is_empty_v<Tag>, "Member tag types are required to be empty");

  public:
    template <typename... Args> AtomMem(Args&&... args) : m_data(std::forward<Args>(args)...) {}

    [[nodiscard]] typename Tag::vector_type const& operator()(Tag) const { return m_data; }
    [[nodiscard]] typename Tag::vector_type& operator()(Tag) { return m_data; }

  private:
    typename Tag::vector_type m_data;
  };

  /**
   * @brief Holds an eigen Array of types defined by the Tag type and provides utilities to fetch
   * the ith column.
   *
   * This is an internal type for use by AtomArray. Each member is tagged with Tag hence, they can
   * be overloaded and selected with the tag type.
   *
   * @tparam Tag An empty type, derived from otf::MemTag.
   */
  template <typename Tag> class EigenArrayAdaptor {
  private:
    static_assert(std::is_empty_v<Tag>, "Member tag types are required to be empty");

    typename Tag::matrix_type m_data;

  public:
    /**
     * @brief Construct a new empty Eigen Array Adaptor object.
     */
    EigenArrayAdaptor() = default;

    /**
     * @brief Construct a new default initialised array with n columns.
     */
    EigenArrayAdaptor(std::size_t n) : m_data(Tag::extent, n) {}

    /**
     * @brief Get the number of elements/columns in the array
     */
    [[nodiscard]] std::size_t size() const noexcept { return m_data.cols(); }

    /**
     * @brief Get the raw {extent by n} underlying array.
     */
    [[nodiscard]] typename Tag::matrix_type const& raw_array(Tag) const noexcept { return m_data; };

    /**
     * @brief Get the raw {extent by n} underlying array.
     */
    [[nodiscard]] typename Tag::matrix_type& raw_array(Tag) noexcept { return m_data; };

    /**
     * @brief Get the ith column of the underlying array, returns a reference to the element in the
     * case that it is a length 1 column.
     */
    [[nodiscard]] decltype(auto) get(Tag, std::size_t i) const {
      if constexpr (Tag::extent == 1) {
        return m_data[i];
      } else {
        return m_data.col(i);
      }
    };

    /**
     * @brief Get the ith column of the underlying array, returns a reference to the element in the
     * case that it is a length 1 column.
     */
    [[nodiscard]] decltype(auto) get(Tag, std::size_t i) {
      if constexpr (Tag::extent == 1) {
        return m_data[i];
      } else {
        return m_data.col(i);
      }
    };
  };

}  // namespace otf::detail