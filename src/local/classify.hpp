#pragma once

#include <vector>

#include "config.hpp"
#include "local/discrete_key.hpp"
#include "local/geometry.hpp"
#include "potentials/neigh_reduce.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"

// Function like object that maps: Universe -> {[discrete_key, ..], [geometry, ...]}.
class Classify {
  public:
    explicit Classify(double r_env) : _r_env(r_env) {}

    void operator()(Supercell const &cell,
                    std::vector<DiscreteKey> &keys,
                    std::vector<Geometry> &geos);

  private:
    struct Index {
        std::size_t idx;
    };

    double _r_env;

    NeighReduce<Index> _reduce;
};

Classify load_classifyer(toml::v2::table const &config);