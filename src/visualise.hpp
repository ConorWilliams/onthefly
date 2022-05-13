#pragma once

#include <optional>

#include "local/classify.hpp"
#include "local/geometry.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"

namespace options {

struct Visualise {
    double r_env = 2.75;  // Radius of environment used to visualise discrete key.

    static Visualise load(toml::v2::table const& config) {
        return {config["visualise"]["r_env"].value_or(Visualise().r_env)};
    }
};

}  // namespace options

class Visualise {
  public:
    explicit Visualise(options::Visualise opt) : _opt(opt), _cl(opt.r_env) {}

    // Colour supercell according to discrete-key of each atom.
    Supercell operator()(Supercell const& cell) {
        // Coloured supercell
        Supercell col = cell;

        _cl(cell, _keys, _geos);

        for (std::size_t i = 0; i < _keys.size(); i++) {
            auto [it, inserted] = _col.try_emplace(_keys[i], _col.size());

            col.activ[i].col.atomic = it->first.sdf[col.activ[i].col];
        }

        return col;
    }

  private:
    options::Visualise _opt;
    Classify _cl;

    std::vector<DiscreteKey> _keys;
    std::vector<Geometry> _geos;

    std::map<DiscreteKey, std::size_t> _col;
};
