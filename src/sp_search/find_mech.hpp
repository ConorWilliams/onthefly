#pragma once

#include <cstdlib>
#include <iostream>
#include <memory>
#include <optional>

#include "config.hpp"
#include "local/environment.hpp"
#include "minimise/LBFGS/lbfgs.hpp"
#include "sp_search/sp_search_base.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"
#include "utility.hpp"
#include "vineyard.hpp"

namespace options {

struct FindMechanisms {
    double r_perturbation;    // (Angstroms) Of perturbation sphere
    double stddev;            // (Angstroms) Of Gaussian perturbations
    double const_pre_factor;  // Ignored if vineyard = true

    double vine_zero_tol = 1e-7;   // Ignored if vineyard = false
    std::size_t consecutive = 10;  // Number of consecutive rediscoveries before finishing
    std::size_t max_search = 50;   // Maximum number of searches
    bool vineyard = false;         // Compute prefactor using vineyard approximation

    Mechanism proto;  // Note: use mech options for Proto-mechanisms

    static FindMechanisms load(toml::v2::table const& config);
};

}  // namespace options

std::vector<ProtoMech> find_mechanisms(options::FindMechanisms const& opt,
                                       Workcell const& init,
                                       std::unique_ptr<PotentialBase>& ff,
                                       std::unique_ptr<SearchBase>& finder);
