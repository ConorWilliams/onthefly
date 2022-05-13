#pragma once

#include <cmath>
#include <cstddef>
#include <iostream>
#include <type_traits>
#include <utility>

#include "config.hpp"
#include "minimise/LBFGS/core.hpp"
#include "minimise/minimiser_base.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

struct MinimiseLBFGS {
    std::size_t n = 10;           // Number of previous steps held in memory
    std::size_t iter_max = 2000;  // Number of translations before exit
    double f2norm = 1e-5;         // Force convergence criterion (eV/Angstroms)
    double proj_tol = 0;          // Trust tolerance
    double max_trust = 0.5;       // Maximum trust radius (Angstroms)
    double min_trust = 0.05;      // Minimum trust radius (Angstroms)
    double grow_trust = 1.5;      // Trust radius expansion rate
    double shrink_trust = 0.5;    // Trust radius contraction rate

    static MinimiseLBFGS load(toml::v2::table const &config);
};

}  // namespace options

// LBFGS optimiser implementation, ignores minimum mode
class MinimiseLBFGS final : public MinimiserBase {
  public:
    explicit MinimiseLBFGS(options::MinimiseLBFGS const &opt) : _opt{opt}, _core{opt.n} {}

    std::unique_ptr<MinimiserBase> clone() const override;

    bool minimise(Supercell &cell,
                  std::unique_ptr<PotentialBase> &ff,
                  VecN<double> const &) override {
        return minimise(cell, ff);
    }

    bool minimise(Supercell &cell, std::unique_ptr<PotentialBase> &ff) override;

  private:
    options::MinimiseLBFGS _opt;
    CoreLBFGS _core;

    VecN<double> _q;
    VecN<double> _gx;
};
