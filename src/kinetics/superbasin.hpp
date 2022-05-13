#pragma once

#include <list>
#include <vector>

#include "config.hpp"
#include "kinetics/basin.hpp"
#include "local/catalogue.hpp"
#include "local/environment.hpp"
#include "local/geometry.hpp"
#include "supercell.hpp"

// Manages a superbasin: a collection of low-barrier-linked basins.
class Superbasin {
  public:
    explicit Superbasin(Basin &&basin) { expand_occupy(std::move(basin)); }

    std::size_t size() const { return _super.size(); }

    // Connect to mechanism 'mech' from basin 'basin' to the currently occupied basin.
    void connect_from(std::size_t basin, std::size_t mech);

    // Expand the superbasin by adding 'basin' to it and setting 'basin' as the currently occupied
    // basin. Returns the previously occupied basin's index.
    std::size_t expand_occupy(Basin &&basin);

    // Search through the basins in the superbasin, if one is found that matches (all active atoms
    // within L2 tolerance 'tol') make it the occupied basin. Returns the previously occupied
    // basin's index.
    std::optional<std::size_t> find_occupy(Supercell const &state, double tol);

    // Deference for current basin, corresponding to the current state of the system
    Basin const &operator*() {
        CHECK(_occupied < size(), "Invalid occupied");
        return _super[_occupied];
    }

    Basin const &operator[](std::size_t i) {
        CHECK(i < size(), "Invalid occupied");
        return _super[i];
    }
    std::size_t occupied() { return _occupied; }

    // Choose a mechanism using the modified mean-rate-method, the mechanism may start from a basin
    // that the system is not currently in, if this is the case the occupied basin is updated and
    // choice.basin_changed is set to true.
    Basin::choice kmc_choice(pcg64 &);

  private:
    std::vector<Basin> _super{};  // Collection of basins
    std::size_t _occupied{};      // Index of current basin in super;
    MatN<double> _prob{};         // Transition probability matrix;

    VecN<double> compute_tau() const;
};
