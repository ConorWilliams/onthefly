#pragma once

#include <cstddef>
#include <memory>

#include "libatom/neighbour/neighbour_list.hpp"
#include "libatom/system/sim_cell.hpp"

namespace otf {

  /**
   * @brief Specifies the virtual-interface for potentials in OLKMC.
   */
  class Potential {
  public:
    // Copy this object
    virtual std::unique_ptr<Potential> clone() const = 0;

    /**
     * @brief Get this potentials cut-off radius that the neighbour lists should be configured for.
     */
    virtual floating rcut() const = 0;

    /**
     * @brief Compute the energy, assumes the neighbour list are ready.
     *
     * Ignores contribution from frozen atoms.
     */
    virtual floating energy(SimCell const &, NeighbourList const &, std::size_t threads) = 0;

    /**
     * @brief Compute gradient, assumes the neighbour list are ready.
     *
     * Force on frozen atoms will be zero.
     */
    virtual void gradient(SimCell &, NeighbourList const &, std::size_t threads) = 0;

    /**
     * @brief Compute hessian matrix, assumes the neighbour list are ready.
     *
     * The resulting hessian will be m by m and only include contributions from the m active atoms.
     */
    virtual void hessian(SimCell const &, NeighbourList const &, std::size_t threads) = 0;

    /**
     * @brief Call parent destructor.
     */
    virtual ~Potential() {}

  protected:
    /**
     * @brief Protected constructor as this is an interface class.
     */
    constexpr Potential() noexcept = default;
  };

}  // namespace otf