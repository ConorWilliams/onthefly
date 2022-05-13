#include "kinetics/basin.hpp"
#include "kinetics/superbasin.hpp"
#include "toml++/toml.h"

namespace options {

struct SuperCache {
    double state_tol;            // (Angstroms) L2 norm between active atoms for states to match.
    double barrier_tol;          // (eV)
    std::size_t cache_size = 8;  // (Superbasins)

    bool dynamic_tol = false;         // If true barrier_tol is dynamically adjusted
    std::size_t max_superbasin_size;  // Maximum number basins in SB before tol_shrink
    double tol_grow;                  // Multiplier by which barrier_tol increases by
    double tol_shrink;                // Multiplier by which barrier_tol decreases by

    Basin basin;

    static SuperCache load(toml::v2::table const &config);
};

}  // namespace options

// Manages a collection of superbasins.
class SuperCache {
  public:
    SuperCache(options::SuperCache opt,
               Supercell const &cell,
               std::vector<Catalogue::pointer> const &env)
        : _opt(opt), _sb({opt.basin, cell, env}) {}

    // Pre condition: _sb's occupied basin's state matches 'cell'.
    // Select a mechanism using appropriate KMC algorithm from active superbasin, if mechanism
    // starts from a different basin cell's active atoms are updated accordingly.
    Basin::choice select_mech(Supercell &cell);

    // Fetch mechanism 'mech' in _sb's active basin.
    Basin::local_mech reconstruct(std::size_t mech) { return (*_sb)[mech]; }

    // Fetch a particular basin in superbasin
    Basin const &operator[](std::size_t i) {
        CHECK(i < _sb.size(), "Invalid occupied");
        return _sb[i];
    }

    // Connect _sb's active basin to the basin 'cell' is currently in via mechanism 'mech', updates
    // active basin to match 'cell'/'env'.
    // Post condition: _sb's occupied basin's state matches 'cell'.
    void connect_via(std::size_t mech,
                     Supercell const &cell,
                     std::vector<Catalogue::pointer> const &env);

    std::size_t size() { return 1 + _cache.size(); }

  private:
    options::SuperCache _opt;

    Superbasin _sb;
    std::list<Superbasin> _cache;

    std::size_t _in_cache_count = 0;

    //////////////////////////////////

    // Caches SB and removes oldest SB.
    void cache(Superbasin &&);
};
