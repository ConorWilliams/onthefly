#include "discrete/discrete_classify.hpp"

#include <vector>

#include "discrete/discrete_supercell.hpp"
#include "local/discrete_key.hpp"
#include "local/geometry.hpp"
#include "potentials/neigh_reduce.hpp"
#include "supercell.hpp"

DiscreteClassify load_discrete_classifyer(toml::v2::table const &config) {
    return DiscreteClassify{options::Discrete::load(config)};
}

// Calculate the nearest neighbours for given lattice point: constant throughout program, though
// occupying atom may change

void DiscreteClassify::initialise(DiscreteSupercell const &cell) {
    _reduce.load(_opt.r_neigh, return_continuous_supercell(cell, true));

    std::size_t i = 0;
    // Label atoms
    for (auto it = _reduce.begin_activ(); it != _reduce.begin_ghost(); ++it) {
        it->idx = i++;
    }

    // Label ghosts
    _reduce.broadcast_ghost_data();

    neigh_list.resize(cell.size());

    for (auto it = _reduce.begin_activ(); it != _reduce.begin_ghost(); ++it) {
        // Reuse memory, reduce allocation
        neigh_list[cell[it->idx].lattice_site].clear();

        _reduce.neigh_reduce(it, [&](auto n) {
            neigh_list[cell[it->idx].lattice_site].emplace_back(cell[n->idx].lattice_site);
        });
    }
}

void DiscreteClassify::operator()(DiscreteSupercell const &cell,
                                  std::vector<DiscreteKey> &keys,
                                  std::vector<Geometry> &geos) {
    geos.resize(cell.vacant.size());
    keys.resize(cell.vacant.size());

    for (std::size_t atom = 0; atom < cell.vacant.size(); atom++) {
        keys[atom].clear();
        geos[atom].clear();

        keys[atom].centre_col = cell.vacant[atom].col;
        keys[atom].sdf[cell.vacant[atom].col] += 1;
        geos[atom].emplace_back(cell._lattice[cell.vacant[atom].lattice_site],
                                cell.vacant[atom].col,
                                cell.occupied() + atom);

        for (std::size_t neigh : neigh_list[cell.vacant[atom].lattice_site]) {
            keys[atom].sdf[cell[cell._map[neigh]].col]++;
            geos[atom].emplace_back(
                cell._lattice[neigh], cell[cell._map[neigh]].col, cell._map[neigh]);
        }

        geos[atom].finalise();
    }
}
