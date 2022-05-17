#include <cstddef>
#include <type_traits>

#include "libatom/asserts.hpp"
#include "libatom/utils.hpp"

namespace otf::detail {

  /**
   * @brief Holds an eigen Array of types defined by the Tag type and provides utilities to fetch
   * the ith column.
   *
   * This is an internal type for use by AtomVector. Each member is tagged with Tag hence, they can
   * be overloaded and selected with the tag type.
   *
   * @tparam Tag, an empty type, derived from otf::Member, defining a scalar_type and extent.
   */
  template <typename Tag> class EigenArrayAdaptor {
  private:
    static_assert(std::is_empty_v<Tag>, "Tag (Member) types are required to be empty");

    Eigen::Array<typename Tag::scalar_type, Tag::extent, Eigen::Dynamic> m_data;

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
     * @brief Get the number of ellements/columns in the array
     */
    std::size_t size() const noexcept { return m_data.cols(); }

    /**
     * @brief Get the raw {extent by n} underlying array.
     */
    auto const &raw_array(Tag) const noexcept { return m_data; };

    /**
     * @brief Get the raw {extent by n} underlying array.
     */
    auto &raw_array(Tag) noexcept { return m_data; };

    /**
     * @brief Get the ith column of the underlying array, returns a reference to the element in the
     * case that it is a length 1 column.
     */
    decltype(auto) get(Tag, std::size_t i) const {
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
    decltype(auto) get(Tag, std::size_t i) {
      if constexpr (Tag::extent == 1) {
        return m_data[i];
      } else {
        return m_data.col(i);
      }
    };
  };

}  // namespace otf::detail