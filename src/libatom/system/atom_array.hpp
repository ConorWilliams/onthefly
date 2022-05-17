#pragma once

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <type_traits>

#include "libatom/asserts.hpp"
#include "libatom/system/detail/eigen_array_adaptor.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief Models an array of "atoms"
   *
   * The default container type used in the libatom; an AtomArray models an array of * "atom" types
   * but decomposes the atom type and stores each member in a separate array. This enables efficient
   * cache use. The members of the "atom" are described through a series of template parameters
   * which should inherit from otf::Member. A selection of canonical members are provided in
   * libatom/system/member.hpp. The members of each atom can be accessed either by the index of the
   * atom or as an eigen array to enable collective operations.
   *
   * Example of use:
   *
   * @code{.cpp}
   *
   * #include "libatom/system/atom_array.hpp"
   * #include "libatom/system/member.hpp"
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
      ((raw_array(Mems{}).conservativeResize(Eigen::NoChange, n)), ...);
    }

    /**
     * @brief Destructive resize, just allocates new memeory
     */
    void destructive_resize(std::size_t n) {
      ((raw_array(Mems{}).conservativeResize(Eigen::NoChange, n)), ...);
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
     * @brief Fetch the entire {extent by size()} underlying array of the Tag member.
     */
    template <typename Tag> [[nodiscard]] auto operator()(Tag) const {
      // Return .leftCols as this is a view type which prevents resizing but allows modification of
      // the values.
      return raw_array(Tag{}).leftCols(size());
    }

    /**
     * @brief Fetch the entire {extent by size()} underlying array of the Tag member.
     */
    template <typename Tag> [[nodiscard]] auto operator()(Tag) {
      // Return .leftCols as this is a view type which prevents resizing but allows modification of
      // the values.
      return raw_array(Tag{}).leftCols(size());
    }

  private:
    static_assert(sizeof...(Mems) > 0, "Need at least one member in an AtomArray");

    using detail::EigenArrayAdaptor<Mems>::raw_array...;
    using detail::EigenArrayAdaptor<Mems>::get...;
  };

}  // namespace otf