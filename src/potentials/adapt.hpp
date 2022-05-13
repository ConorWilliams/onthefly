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
#include "utility.hpp"

// Compute gradient of norm_sq(grad(u))
template <typename Potential> class Adapt final : public PotentialBase {
  public:
    Adapt(Adapt const &) = default;

    template <typename... Args> Adapt(Args &&...args) : _nabla(std::forward<Args>(args)...) {}

    std::unique_ptr<PotentialBase> clone() const override {
        return std::make_unique<Adapt<Potential>>(*this);
    }

    // Get this force fields cut-off radius
    double rcut() const override { return _nabla.rcut(); }

    std::unordered_map<std::string, std::uint16_t> const &species_map() const override {
        return _nabla.species_map();
    }

    // Compute energy
    double energy(Supercell const &) override { ALWAYS_CHECK(false, "Grad not supported"); }

    // Compute gradient of norm_sq(grad(u))
    void gradient(Supercell const &cell, VecN<double> &out) override {
        _x = cell;
        _nabla.gradient(_x, _g);

        double delta = 0.01 / norm(_g);

        _x.activ.view() += delta * _g;
        _nabla.gradient(_x, out);

        out = (out - _g) * (1 / delta);
    }

    // Call parent destructor
    virtual ~Adapt() {}

  private:
    Potential _nabla;
    Supercell _x;
    VecN<double> _g;
};