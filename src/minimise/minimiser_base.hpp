#pragma once

#include <memory>
#include <random>

#include "config.hpp"
#include "pcg_random.hpp"
#include "potentials/potential_base.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

// Produce a random vector on the n dimensional sphere
template <typename T> void rand_normal(T &x) {
    //
    static thread_local pcg64 rng(pcg_extras::seed_seq_from<std::random_device>{});

    std::normal_distribution<> gauss_dist(0, 1);

    for (auto &&elem : x) {
        elem = gauss_dist(rng);
    }

    x *= 1 / norm(x);
}

// Here we define the virtual-interface for minimisers in OLKMC,
class MinimiserBase {
  public:
    // Return a copy of this
    virtual std::unique_ptr<MinimiserBase> clone() const = 0;

    // Move supercell to local minimum of potential, vector is an approximation of min-mode
    virtual bool minimise(Supercell &, std::unique_ptr<PotentialBase> &, VecN<double> const &) = 0;

    // Move supercell to local minimum of potential
    virtual bool minimise(Supercell &, std::unique_ptr<PotentialBase> &) = 0;

    // Call parent destructor
    virtual ~MinimiserBase() = default;

  protected:
    // Protected constructor as this is an interface class
    constexpr MinimiserBase() noexcept = default;
};

// Customisation point, dynamically select minimiser
std::unique_ptr<MinimiserBase> load_minimiser(toml::v2::table const &config);