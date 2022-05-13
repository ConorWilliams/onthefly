#pragma once

#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "config.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"

// Here we define the virtual-interface for potentials in OLKMC,
class PotentialBase {
  public:
    // Copy this object
    virtual std::unique_ptr<PotentialBase> clone() const = 0;

    // Get this potentials cut-off radius
    virtual double rcut() const = 0;

    // Maps species string to index this potential uses to internally represent that species
    virtual std::unordered_map<std::string, std::uint16_t> const &species_map() const = 0;

    // Compute energy
    virtual double energy(Supercell const &cell) = 0;

    // Compute gradient
    virtual void gradient(Supercell const &cell, VecN<double> &out) = 0;

    // Compute mass-weighted hessian
    virtual void hessian(Supercell const &, MatN<double> &) {
        throw std::runtime_error("Unimplemented");
    };

    // Call parent destructor
    virtual ~PotentialBase() {}

  protected:
    // Protected constructor as this is an interface class
    constexpr PotentialBase() noexcept = default;
};

// Customisation point, dynamically select potential
std::unique_ptr<PotentialBase> load_potential(toml::v2::table const &config);
