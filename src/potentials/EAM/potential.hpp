#pragma once

#include <cstddef>
#include <memory>
#include <unordered_map>

#include "config.hpp"
#include "potentials/EAM/data.hpp"
#include "potentials/neigh_reduce.hpp"
#include "potentials/potential_base.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"

// Force functor for any tabulated EAM potential, holds a shared pointer to EAM
// data therefore, use copy constructor when passing to threads
class PotentialEAM final : public PotentialBase {
  public:
    explicit PotentialEAM(TabEAM &&data);

    std::unique_ptr<PotentialBase> clone() const override;

    // Get this force fields cut-off radius
    double rcut() const override { return _data->rcut(); }

    std::unordered_map<std::string, std::uint16_t> const &species_map() const override {
        return _data->species_map();
    }

    // Compute energy
    double energy(Supercell const &cell) override;

    // Compute gradient
    void gradient(Supercell const &cell, VecN<double> &out) override;

    // Compute gradient
    void hessian(Supercell const &cell, MatN<double> &out) override;

  private:
    struct Empty {};

    struct Grad {
        double fp_rho = 0;
    };
    struct Hess {
        double rho = 0;
        Vec3<double> mu = Vec3<double>::Zero();
        std::size_t idx = 0;
    };

    std::shared_ptr<TabEAM const> _data;

    NeighReduce<Empty> _reduce_energy;
    NeighReduce<Grad> _reduce_grad;
    NeighReduce<Hess> _reduce_hess;
};
