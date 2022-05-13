#include "local/classify.hpp"

#include <vector>

#include "local/discrete_key.hpp"
#include "local/geometry.hpp"
#include "potentials/neigh_reduce.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"

Classify load_classifyer(toml::v2::table const &config) {
    return Classify{fetch<double>(config, "catalogue", "r_env")};
}

// Maps: Universe -> {[discrete_key, ..], [geometry, ...]}.
void Classify::operator()(Supercell const &cell,
                          std::vector<DiscreteKey> &keys,
                          std::vector<Geometry> &geos) {
    //
    _reduce.load(_r_env, cell);

    std::size_t i = 0;
    // Label atoms
    for (auto it = _reduce.begin_activ(); it != _reduce.begin_bound(); ++it) {
        it->idx = i++;
    }
    // Label ghosts
    _reduce.broadcast_ghost_data();

    keys.resize(cell.activ.size());  // Usually a no-op
    geos.resize(cell.activ.size());  // Usually a no-op

    for (auto it = _reduce.begin_activ(); it != _reduce.begin_bound(); ++it) {
        // Resuse memory, reduce allocation
        keys[it->idx].clear();
        geos[it->idx].clear();

        // As reduce does not include it/central-atom
        keys[it->idx].centre_col = it->col;
        keys[it->idx].sdf[it->col] += 1;
        geos[it->idx].emplace_back(it->vec, it->col, it->idx);

        // Add neighbours
        _reduce.neigh_reduce(it, [&](auto n, double, Eigen::Array3d const &) {
            keys[it->idx].sdf[n->col]++;
            geos[it->idx].emplace_back(n->vec, n->col, n->idx);
        });

        geos[it->idx].finalise();
    }
}
