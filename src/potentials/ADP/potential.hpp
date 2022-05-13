#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "potentials/ADP/data.hpp"
#include "potentials/neigh_reduce.hpp"
#include "potentials/potential_base.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"

// Force functor for any tabulated ADP potential, holds a shared pointer to ADP
// data therefore, use copy constructor when passing to threads
class PotentialADP final : public PotentialBase {
  public:
    explicit PotentialADP(TabADP &&data);

    std::unique_ptr<PotentialBase> clone() const override;

    std::unordered_map<std::string, std::uint16_t> const &species_map() const override {
        return _data->species_map();
    }

    double rcut() const override { return _data->rcut(); }

    // Compute energy
    double energy(Supercell const &cell) override;

    // Compute force
    void gradient(Supercell const &cell, VecN<double> &out) override;

  private:
    struct Density {
        double density = 0;
        Vec3<double> dipole_density = Vec3<double>::Zero();
        Mat3<double> quadropole_density = Mat3<double>::Zero();
    };

    std::shared_ptr<TabADP const> _data;

    NeighReduce<Density> _reduce;
};
