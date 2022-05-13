#pragma once

#include <limits>

#include "config.hpp"
#include "local/catalogue.hpp"
#include "local/environment.hpp"
#include "local/geometry.hpp"
#include "pcg_random.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"

namespace options {

struct Basin {
    double temp;  // k
    double max_barrier = std::numeric_limits<double>::max();
    // eV Any mechanisms with an activation energy higher than this will be ignored

    static Basin load(toml::v2::table const &config);
};

}  // namespace options

// Represents a basin of the potential energy of the entire system. Stores a reference image of the
// system (active atoms only) and a list of all mechanisms accessible from the basin. Implements
// standard kmc algorithm on this list. Stored mechanisms are pointers into catalogue and can be
// reconstructed onto a supercell.
class Basin {
  public:
    struct choice {
        bool basin_changed;  // True if mech starts from different basin.
        std::size_t mech;    // i'th mech in current basin.
        double delta_t;      // Time increment if mech is carried out.
        std::size_t basin;   // Current basin in superbasin
    };

    // Represents generalised-mechanism acting on a SPECIFIC atom.
    class local_mech {
      public:
        double rate;            // (Hz)
        double barrier;         // Maximum of forward/reverse barriers.
        bool exit_mech = true;  // If true it is known *this links: inside SB -> outside SB.

        local_mech(double rate,
                   double barrier,
                   std::size_t atom_idx,
                   std::size_t mech_off,
                   Catalogue::pointer env)
            : rate(rate), barrier(barrier), _env(env), _atom_idx(atom_idx), _mech_off(mech_off) {}

        // Reconstruct mechanism pointed to by *this onto supercell, returns ref to mech
        // reconstructed.
        Mechanism const &onto(Supercell &, std::vector<Geometry> &) const;

        void refine(std::vector<Geometry> &geo) const;

      private:
        Catalogue::pointer _env;
        std::size_t _atom_idx;
        std::size_t _mech_off;
    };

    Basin(options::Basin const &, Supercell const &, std::vector<Catalogue::pointer> const &);

    // Uses standard n-fold-way kmc algorithm to select mechanism inside basin. choice.basin_changed
    // is always false,
    choice kmc_choice(pcg64 &, std::size_t basin) const;

    local_mech &operator[](std::size_t i) {
        CHECK(i < size(), "Invalid mech");
        return _mechs[i];
    }

    local_mech const &operator[](std::size_t i) const {
        CHECK(i < size(), "Invalid mech");
        return _mechs[i];
    }

    std::size_t size() const { return _mechs.size(); }

    VecN<double> const &state() const { return _state; }
    double rate_sum() const { return _rate_sum; }

  public:
    bool connected = false;  // True if any of the mechs have exit_mech = false

  private:
    VecN<double> _state;
    std::vector<local_mech> _mechs;
    double _rate_sum = 0;
};
