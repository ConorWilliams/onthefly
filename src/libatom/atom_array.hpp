#pragma once

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <type_traits>

#include "libatom/detail/eigen_array_adaptor.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief A base type to derive from for defining members of an atom for use in AtomArrays.
   *
   * @tparam Scalar The type of the atoms member.
   * @tparam Extent How many elements of the Scalar type are in the member (vector dimension).
   */
  template <typename Scalar, std::size_t Extent = 1> struct AtomArrayMem {
    //
    static_assert(std::is_default_constructible_v<Scalar>);

    /** @brief This member represents a vector of elements of scalar_type. */
    using scalar_type = Scalar;
    /** @brief A vector of scalar_type length extent. */
    using vector_type = std::conditional_t<Extent == 1, Scalar, Eigen::Array<Scalar, 1, Extent>>;
    /** @brief The underlying Eigen type used to store a dynamic number of these members. */
    using matrix_type = Eigen::Array<Scalar, Extent, Eigen::Dynamic>;

    static constexpr std::size_t extent = Extent;
  };

  /**
   * @brief Models an array of "atoms"
   *
   * The default container type used in the libatom; an AtomArray models an array of * "atom" types
   * but decomposes the atom type and stores each member in a separate array. This enables efficient
   * cache use. The members of the "atom" are described through a series of template parameters
   * which should inherit from otf::AtomArrayMem. A selection of canonical members are provided in
   * libatom/member.hpp. The members of each atom can be accessed either by the index of the
   * atom or as an eigen array to enable collective operations.
   *
   * Example of use:
   *
   * @code{.cpp}
   *
   * #include "libatom/atom_array.hpp"
   *
   * using namespace otf;
   *
   * AtomArray<Pos, AtomicNum> atoms{10}; // Initialise an array of 10 atoms.
   *
   * // Add a hydrogen atom at the origin.
   *
   * atoms(Pos{}, 0) =  Vec3<double>{0, 0, 0};
   * atoms(AtomicNum{}, 0) = 1;
   *
   * Vec3 xyz = atoms(Pos{}, 0); // Get the position of the zeroth atom.
   *
   * std::size_t n = = atoms(AtomicNum{}, 0); // Get the atomic number of the zeroth atom.
   *
   * atoms(Pos{}) += 1; // Add 1 to each of every atoms coordinates.
   *
   * atoms(AtomicNum{}) = 6; // Set all the atoms to carbon atoms.
   *
   * @endcode
   */
  template <typename... Mems> class AtomArray : private detail::EigenArrayAdaptor<Mems>... {
  public:
    /**
     * @brief Construct a new empty AtomArray.
     */
    AtomArray() = default;

    /**
     * @brief Construct an atom vector of n atoms with all default initialised members.
     */
    explicit AtomArray(std::size_t n) : detail::EigenArrayAdaptor<Mems>(n)... {}

    /**
     * @brief Fetch the used size of each array
     */
    using detail::EigenArrayAdaptor<typename First<Mems...>::type>::size;

    /**
     * @brief Constructive resize, allocates+copies unless size is the same
     */
    void resize(std::size_t n) {
      ((void)(raw_array(Mems{}).conservativeResize(Eigen::NoChange, n)), ...);
    }

    /**
     * @brief Destructive resize, just allocates new memeory
     */
    void destructive_resize(std::size_t n) {
      ((void)(raw_array(Mems{}).conservativeResize(Eigen::NoChange, n)), ...);
    }

    /**
     * @brief Fetch the ith element of the Tag member.
     *
     * @return decltype(auto) Either an Eigen view into a vector or a reference to the element if
     * 1D.
     */
    template <typename Tag> [[nodiscard]] decltype(auto) operator()(Tag, std::size_t i) const {
      return get(Tag{}, i);
    }

    /**
     * @brief Fetch the ith element of the Tag member.
     *
     * @return decltype(auto) Either an Eigen view into a vector or a reference to the element if
     * 1D.
     */
    template <typename Tag> [[nodiscard]] decltype(auto) operator()(Tag, std::size_t i) {
      return get(Tag{}, i);
    }

    /**
     * @brief Fetch a view of the entire {extent by size()} underlying array of the Tag member.
     */
    template <typename Tag> [[nodiscard]] auto operator()(Tag) const {
      return raw_array(Tag{}).leftCols(size());
    }

    /**
     * @brief Fetch a view of the entire {extent by size()} underlying array of the Tag member.
     */
    template <typename Tag> [[nodiscard]] auto operator()(Tag) {
      return raw_array(Tag{}).leftCols(size());
    }

  private:
    static_assert(sizeof...(Mems) > 0, "Need at least one member in an AtomArray");

    using detail::EigenArrayAdaptor<Mems>::raw_array...;
    using detail::EigenArrayAdaptor<Mems>::get...;
  };

  /**
   * @brief A collection of default member types for use in AtomArray's
   */
  namespace members {

    /**
     * @brief Tag type for position (xyz).
     */
    struct Position : AtomArrayMem<floating, spatial_dims> {};

    /**
     * @brief Tag type for dimer axis (xyz).
     */
    struct Axis : AtomArrayMem<floating, spatial_dims>{};

    /**
     * @brief Tag type for gradiant of the potential.
     */
    struct Gradient : AtomArrayMem<floating, spatial_dims> {};

    /**
     * @brief Tag type for atomic number.
     */
    struct Velocity : AtomArrayMem<floating, spatial_dims> {};

    /**
     * @brief Tag type for atomic number.
     */
    struct AtomicNum : AtomArrayMem<std::size_t, 1> {};

    /**
     * @brief Tag type for atomic number.
     */
    struct Mass : AtomArrayMem<std::size_t, 1> {};

    /**
     * @brief Tag type for index.
     */
    struct Index : AtomArrayMem<std::size_t, 1> {};

    /**
     * @brief Tag type for index.
     */
    struct Symbol : AtomArrayMem<std::string_view, 1> {};

    /**
     * @brief Tag type for frozen atoms.
     */
    struct Frozen : AtomArrayMem<bool, 1> {};

  }  // namespace members

  // No online namespace as m.css doens't like them :()
  using namespace members;

}  // namespace otf