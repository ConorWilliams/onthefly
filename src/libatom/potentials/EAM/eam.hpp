#pragma once

#include <cstddef>
#include <memory>

#include "libatom/asserts.hpp"
#include "libatom/neighbour/neighbour_list.hpp"
#include "libatom/potentials/EAM/data.hpp"
#include "libatom/potentials/potential.hpp"
#include "libatom/system/atom_array.hpp"
#include "libatom/system/member.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf {

  class EAM : public Potential {
  public:
    EAM(std::shared_ptr<DataEAM const> data) : m_data(std::move(data)) {}

    /**
     * @brief Copies this potentials cut-off radius
     */
    std::unique_ptr<Potential> clone() const override { return std::make_unique<EAM>(*this); }

    /**
     * @brief Get this potentials cut-off radius that the neighbour lists should be configured for.
     */
    floating rcut() const override { return m_data->rcut(); }

    /**
     * @brief Compute the energy, assumes the neighbour list are ready.
     *
     * Ignores contribution from frozen atoms.
     */
    floating energy(SimCell const &, NeighbourList const &, std::size_t num_threads) override;

    /**
     * @brief Compute gradient, assumes the neighbour list are ready.
     *
     * Force on frozen atoms will be zero.
     */
    void gradient(SimCell &, NeighbourList const &, std::size_t num_threads) override;

    /**
     * @brief Compute hessian matrix, assumes the neighbour list are ready.
     *
     * The resulting hessian will be m by m and only include contributions from the m active atoms.
     */
    void hessian(SimCell const &, NeighbourList const &, std::size_t num_threads) override;

  private:
    std::shared_ptr<DataEAM const> m_data;

    struct Fprime : Member<floating, 1> {};
    struct Rho : Member<floating, 1> {};
    struct Mu : Member<floating, 3> {};

    AtomArray<Fprime, Rho, Mu> m_aux;
  };

}  // namespace otf