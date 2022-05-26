#pragma once

#include <cstddef>
#include <memory>

#include "libatom/asserts.hpp"
#include "libatom/atom_array.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/potentials/EAM/data.hpp"
#include "libatom/potentials/base.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf::potentials {

  class EAM : public Base {
  public:
    EAM(std::shared_ptr<DataEAM const> data) : m_data(std::move(data)) {}

    /**
     * @brief Copies this potentials cut-off radius
     */
    std::unique_ptr<Base> clone() const override { return std::make_unique<EAM>(*this); }

    /**
     * @brief Get this potentials cut-off radius that the neighbour lists should be configured for.
     */
    floating rcut() const override { return m_data->rcut(); }

    /**
     * @brief Compute the energy, assumes the neighbour list are ready.
     *
     * Ignores contribution from frozen atoms.
     */
    floating energy(SimCell const &, neighbour::List const &, std::size_t num_threads) override;

    /**
     * @brief Compute gradient, assumes the neighbour list are ready.
     *
     * Force on frozen atoms will be zero.
     */
    void gradient(SimCell &, neighbour::List const &, std::size_t num_threads) override;

    /**
     * @brief Compute hessian matrix, assumes the neighbour list are ready.
     *
     * The resulting hessian will be m by m and only include contributions from the m active atoms.
     */
    void hessian(SimCell const &, neighbour::List const &, std::size_t num_threads) override;

  private:
    std::shared_ptr<DataEAM const> m_data;

    struct Fprime : AtomArrayMem<floating, 1> {};
    struct Rho : AtomArrayMem<floating, 1> {};
    struct Mu : AtomArrayMem<floating, 3> {};

    AtomArray<Fprime, Rho, Mu> m_aux;
  };

}  // namespace otf::potentials