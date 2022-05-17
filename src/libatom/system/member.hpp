

#include <cstddef>
#include <type_traits>

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

}  // namespace otf