#pragma once

#include <cstddef>
#include <memory>

#include "libatom/asserts.hpp"
#include "libatom/neighbour/neighbour_list.hpp"
#include "libatom/potentials/EAM/data.hpp"
#include "libatom/system/atom_array.hpp"
#include "libatom/system/member.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf {

  class EAM {
  public:
    EAM(std::shared_ptr<DataEAM const> data) : m_data(std::move(data)) {}

    /**
     * @brief Get this potentials cut-off radius
     */
    floating rcut() const { return m_data->rcut(); }

    /**
     * @brief Compute the energy, assumes the neighbour list are ready
     */
    floating energy(SimCell const &, NeighbourList const &, std::size_t num_threads = 1) const;

    // Compute gradient, force on frozen atoms must be zero
    void gradient(SimCell &, NeighbourList const &, std::size_t num_threads = 1);

    // Compute gradient, force on frozen atoms must be zero
    void hessian(SimCell const &, NeighbourList const &, std::size_t num_threads = 1) const;

  private:
    std::shared_ptr<DataEAM const> m_data;

    struct Fprime : Member<floating, 1> {};
    struct Rho : Member<floating, 1> {};
    struct Mu : Member<floating, 3> {};

    AtomArray<Fprime, Rho, Mu> m_aux;
  };

}  // namespace otf