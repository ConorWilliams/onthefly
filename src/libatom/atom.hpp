#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <type_traits>

#include "libatom/asserts.hpp"
#include "libatom/detail/atom.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief A base type to derive from for defining members of an Atom type.
   *
   * @tparam Scalar The type of the atoms member.
   * @tparam Extent How many elements of the Scalar type are in the member (vector dimension).
   */
  template <typename Scalar, std::size_t Extent = 1> struct MemTag {
    //
    static_assert(std::is_default_constructible_v<Scalar>);

    /** @brief This member represents a vector of elements of scalar_type. */
    using scalar_type = Scalar;
    /** @brief A vector of scalar_type length extent. */
    using vector_type = std::conditional_t<Extent == 1, Scalar, Eigen::Array<Scalar, Extent, 1>>;
    /** @brief The underlying Eigen type used to store a dynamic number of these members. */
    using matrix_type = Eigen::Array<Scalar, Extent, Eigen::Dynamic>;

    static constexpr std::size_t extent = Extent;
  };

  /**
   * @brief The otf representation of an atom.
   *
   * @tparam Mems a series of empty types, derived from otf::MemTag, to describe each member.
   */
  template <typename... Mems> struct Atom : detail::AtomMem<Mems>... {
    // Expose tagged dispatch.
    using detail::AtomMem<Mems>::operator()...;

    // /**
    //  * @brief Construct a new Atom object, explicitly defaulted.
    //  */
    // Atom(Atom const&) = default;

    // /**
    //  * @brief Construct a new Atom object, explicitly defaulted.
    //  */
    // Atom(Atom&&) = default;

    /**
     * @brief Construct a new Atom object, forwards each argument to a member.
     */
    template <typename... Args, std::enable_if_t<sizeof...(Mems) == sizeof...(Args), int> = 0>
    Atom(Args&&... args) : detail::AtomMem<Mems>(std::forward<Args>(args))... {}
  };

  /**
   * @brief A container that models a std::vector of Atom types.
   *
   * The members of the "atom" are described through a series of template parameters which should
   * inherit from otf::MemTag. A selection of canonical members are provided in the namespace
   * builtin_members.
   *
   * Example of use:
   *
   * @code{.cpp}
   *
   * #include "libatom/atom.hpp"
   *
   * using namespace otf;
   *
   * AtomVector<Position> atoms; // Initialise a vector of 0 atoms.
   *
   * atoms.emplace_back({1,2,3});
   *
   *
   * @endcode
   */
  template <typename... Mems> class AtomVector : private std::vector<Atom<Mems...>> {
  private:
    static_assert(sizeof...(Mems) > 0, "Need at least one member in an AtomVector");

    using Vector = std::vector<Atom<Mems...>>;

  public:
    // Expose subset of underlying vector API
    using Vector::begin;
    using Vector::clear;
    // using Vector::emplace_back;
    using Vector::end;
    using Vector::push_back;
    using Vector::size;
    using Vector::Vector;

    /**
     * @brief Bounds checked version of operator [].
     */
    decltype(auto) operator[](std::size_t i) {
      ASSERT(i < size(), "Out of bounds");
      return Vector::operator[](i);
    }

    /**
     * @brief Bounds checked version of operator [] const.
     */
    decltype(auto) operator[](std::size_t i) const {
      ASSERT(i < size(), "Out of bounds");
      return Vector::operator[](i);
    }

    /**
     * @brief Provides an emplace_back with explicit vector_types.
     */
    decltype(auto) emplace_back(typename Mems::vector_type const&... args) {
      return Vector::emplace_back(args...);
    }
  };

  /**
   * @brief A container that stores a Atom types decomposed by member.
   *
   * The default container type used in the libatom; an AtomArray models an array of "atom" types
   * but decomposes the atom type and stores each member in a separate array. This enables efficient
   * cache use. The members of the "atom" are described through a series of template parameters
   * which should inherit from otf::MemTag. A selection of canonical members are provided in the
   * namespace builtin_members. The members of each atom can be accessed either by the index of the
   * atom or as an Eigen array to enable collective operations.
   *
   * Example of use:
   *
   * @code{.cpp}
   *
   * #include "libatom/atom.hpp"
   *
   * using namespace otf;
   *
   * AtomArray<Position, AtomicNum> atoms{10}; // Initialise an array of 10 atoms.
   *
   * // Add a hydrogen atom at the origin.
   *
   * atoms(Position{}, 0) =  Vec3<double>{0, 0, 0};
   * atoms(AtomicNum{}, 0) = 1;
   *
   * Vec3 xyz = atoms(Position{}, 0); // Get the position of the zeroth atom.
   *
   * std::size_t n = = atoms(AtomicNum{}, 0); // Get the atomic number of the zeroth atom.
   *
   * atoms(Position{}) += 1; // Add 1 to each of every atoms coordinates.
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
  namespace builtin_members {

    /**
     * @brief Tag type for position (xyz).
     */
    struct Position : MemTag<floating, spatial_dims> {};

    /**
     * @brief Tag type for dimer axis (xyz).
     */
    struct Axis : MemTag<floating, spatial_dims> {};

    /**
     * @brief Tag type for gradient of the potential.
     */
    struct Gradient : MemTag<floating, spatial_dims> {};

    /**
     * @brief Tag type for velocity.
     */
    struct Velocity : MemTag<floating, spatial_dims> {};

    /**
     * @brief Tag type for atomic number.
     */
    struct AtomicNum : MemTag<std::size_t, 1> {};

    /**
     * @brief Tag type for atomic mass.
     */
    struct Mass : MemTag<std::size_t, 1> {};

    /**
     * @brief Tag type for index.
     */
    struct Index : MemTag<std::size_t, 1> {};

    /**
     * @brief Tag type for atomic symbol.
     */
    struct Symbol : MemTag<std::string_view, 1> {};

    /**
     * @brief Tag type for frozen atoms.
     */
    struct Frozen : MemTag<bool, 1> {};

    /**
     * @brief Tag type for atom colour (generalisation of atomic number)
     */
    struct Colour : MemTag<std::size_t, 1> {};

  }  // namespace builtin_members

  // No online namespace as m.css doesn't like them :(
  using namespace builtin_members;

}  // namespace otf