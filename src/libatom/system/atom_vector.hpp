
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <type_traits>

#include "libatom/asserts.hpp"
#include "libatom/system/detail/eigen_array_adaptor.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief A base type to derive from for defining members of an atom for use in AtomVectors.
   *
   * @tparam Scalar The type of the atoms member
   * @tparam Extent How many elements of the Scalar type are in the member (vector dimension)
   */
  template <typename Scalar, std::size_t Extent = 1> struct Member {
    //
    using scalar_type = Scalar;
    using vector_type = std::conditional_t<Extent == 1, Scalar, Eigen::Array<Scalar, 1, Extent>>;

    static constexpr std::size_t extent = Extent;

    static_assert(std::is_default_constructible_v<Scalar>);
  };

  /**
   * @brief Tag type for position (xyz)
   */
  struct Pos : Member<floating, spatial_dims> {};
  /**
   * @brief Tag type for atomic number
   */
  struct AtomicNum : Member<std::size_t, 1> {};

  /**
   * @brief A dynamic array of "atoms"
   *
   * The default adaptive container type used in the libatom; an AtomVector models a vector of
   * "atom" types but decomposes the atom type and stores each member in a separate vector. This
   * enables efficient cache use. The members of the "atom" are described through a series of
   * template parameters which should inherit from otf::Member. A selection of canonical members are
   * provided in libatom/system/member.hpp. The members of each atom can be accessed either by
   * the index of the atom or as an eigen array to enable collective operations.
   *
   * Example of use:
   *
   * @code{.cpp}
   *
   * #include "libatom/system/atomvector.hpp"
   *
   * using namespace otf;
   *
   * AtomVector<Pos, AtomicNum> atoms;
   *
   * atoms.push_back({0,0,0}, 1) // Add a hydrogen atom to the origin
   *
   * Vec3 xyz = atoms(Pos{}, 0) // Get the position of the zeroth atom
   *
   * std::size_t n = = atoms(AtomicNum{}, 0) // Get the atomic number of the zeroth atom
   *
   * atoms(Pos{}) += 1. // Add 1 to each of every atoms coordinates
   *
   * @endcode
   */
  template <typename... Mems> class AtomVector : private detail::EigenArrayAdaptor<Mems>... {
  public:
    /**
     * @brief Construct a new empty AtomVector.
     */
    AtomVector() = default;

    /**
     * @brief Construct an atom vector of n atoms with all default initialised members.
     */
    explicit AtomVector(std::size_t n) : detail::EigenArrayAdaptor<Mems>(n)... {}

    AtomVector(AtomVector &&) = default;

    /**
     * @brief Fetch the used size of each array
     */
    [[nodiscard]] std::size_t size() const noexcept { return m_size; }

    /**
     * @brief Shrink the arrays. Does not reallocate.
     */
    void shrink(std::size_t n) {
      VERIFY(n <= size(), "Cannot 'shrink' larger!");
      m_size = n;
    }
    /**
     * @brief Fetch the ith element of the Tag member.
     *
     * @return decltype(auto) Either an Eigen view into a vector or a reference to the element if
     * 1D.
     */
    template <typename Tag> [[nodiscard]] decltype(auto) operator()(Tag, std::size_t i) const {
      ASSERT(i < size(), "Out of bounds");
      return get(Tag{}, i);
    }

    /**
     * @brief Fetch the ith element of the Tag member.
     *
     * @return decltype(auto) Either an Eigen view into a vector or a reference to the element if
     * 1D.
     */
    template <typename Tag> [[nodiscard]] decltype(auto) operator()(Tag, std::size_t i) {
      ASSERT(i < size(), "Out of bounds");
      return get(Tag{}, i);
    }

    /**
     * @brief Fetch the entire {extent by size()} underlying array of the Tag member.
     */
    template <typename Tag> [[nodiscard]] auto operator()(Tag) const {
      return raw_array(Tag{}).leftCols(size());
    }

    /**
     * @brief Fetch the entire {extent by size()} underlying array of the Tag member.
     */
    template <typename Tag> [[nodiscard]] auto operator()(Tag) {
      return raw_array(Tag{}).leftCols(size());
    }

    /**
     * @brief Add a new element to the back of each array.
     */
    void push_back(typename Mems::vector_type const &...args) { emplace_back(args...); }

    /**
     * @brief Emplac a  new element to the back of each array.
     */
    template <typename... Args> void emplace_back(Args &&...args) {
      //
      static_assert(sizeof...(Args) == sizeof...(Mems), "Too few args to emplace_back");

      ASSERT(size() <= capacity(), "Pre check invariant");

      if (size() == capacity()) {
        // Need to grow arrays
        std::size_t new_cap = (3 * (capacity() + 1)) / 2;

        ((raw_array(Mems{}).conservativeResize(Eigen::NoChange, new_cap)), ...);
      }

      ASSERT(size() < capacity(), "Not enough space");

      ((get(Mems{}, size()) = std::forward<Args>(args)), ...);

      m_size += 1;

      ASSERT(size() <= capacity(), "Pre check invariant");
    }

  private:
    static_assert(sizeof...(Mems) > 0, "Need at least one member in an AtomVector");

    std::size_t m_size = 0;

    using detail::EigenArrayAdaptor<Mems>::raw_array...;
    using detail::EigenArrayAdaptor<Mems>::get...;

    /**
     * @brief Utility to extract first type in a parameter pack.
     */
    template <typename T, typename...> struct First { using type = T; };

    /**
     * @brief Size of the underlying arrays.
     */
    [[nodiscard]] std::size_t capacity() const {
      return detail::EigenArrayAdaptor<typename First<Mems...>::type>::size();
    }
  };

}  // namespace otf