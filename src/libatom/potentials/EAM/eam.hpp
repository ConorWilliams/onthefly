#pragma once

#include <cstddef>
#include <memory>

#include "libatom/atom.hpp"
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
    floating energy(SimCell const &, neighbour::List &, std::size_t num_thread) override;

    /**
     * @brief Compute gradient, assumes the neighbour list are ready.
     *
     * Force on frozen atoms will be zero.
     */
    std::optional<floating> gradient(SimCell &, neighbour::List &, std::size_t num_thread) override;

    /**
     * @brief Compute hessian matrix, assumes the neighbour list are ready.
     *
     * The resulting hessian will be m by m and only include contributions from the m active atoms.
     */
    void hessian(SimCell &, neighbour::List &, std::size_t num_thread) override;

  private:
    std::shared_ptr<DataEAM const> m_data;

    struct Fprime : MemTag<floating, 1> {};
    struct Rho : MemTag<floating, 1> {};
    struct Mu : MemTag<floating, spatial_dims> {};
    struct Hidx : MemTag<std::size_t, 1> {};

    AtomArray<Fprime, Rho, Mu, Hidx> m_aux;
  };

}  // namespace otf::potentials