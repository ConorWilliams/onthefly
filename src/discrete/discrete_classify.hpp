#pragma once

#include "config.hpp"
#include "discrete/discrete_supercell.hpp"
#include "local/discrete_key.hpp"
#include "local/geometry.hpp"
#include "potentials/neigh_reduce.hpp"
#include "supercell.hpp"

class DiscreteClassify {
  public:
    explicit DiscreteClassify(options::Discrete opt) : _opt(opt) {}

    void initialise(DiscreteSupercell const &cell);
    void operator()(DiscreteSupercell const &cell,
                    std::vector<DiscreteKey> &keys,
                    std::vector<Geometry> &geos);

    std::vector<std::vector<std::size_t>> neigh_list{};

  private:
    struct Index {
        std::size_t idx;
    };

    NeighReduce<Index> _reduce;

    options::Discrete _opt;
};

DiscreteClassify load_discrete_classifyer(toml::v2::table const &config);
