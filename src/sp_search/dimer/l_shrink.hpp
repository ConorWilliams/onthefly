#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>

#include "config.hpp"
#include "minimise/LBFGS/core.hpp"
#include "minimise/minimiser_base.hpp"
#include "potentials/potential_base.hpp"
#include "sp_search/dimer/dimer.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

struct LShrink {
    std::size_t iter_max_tran = 1000;  // Number of translations before exit
    std::size_t BB_method = 1;         // Use BB1 vs BB2 https://doi.org/10.1137/19M1253356

    double delta_t;  // Shrinking rate
    double l0;       // (Angstroms) Dimer length
    double l_min;    // (Angstrom) Min dimer length
    double s_max;    // (Angstrom) max step size
    double f2norm;   // Force convergence criterion (eV/Angstroms)

    double nudge;      // (Angstroms) Distance perturbed along min-mode at SP
    double basin_tol;  // (Angstroms) L2 between active atoms to be considered distinct basins

    static LShrink load(toml::v2::table const &config);
};

}  // namespace options

// Class responsible for moving a dimer (cell+axis) to a near-by saddle-point
class LShrink {
  public:
    LShrink(options::LShrink const &opt) : _opt(opt), _core{5} {}

    bool find_sp(Supercell &cell, VecN<double> &ax, std::unique_ptr<PotentialBase> &ff);

    options::LShrink const &get_opt() const { return _opt; }

  private:
    options::LShrink _opt;

    CoreLBFGS _core;

    VecN<double> _x;
    VecN<double> _gx;

    VecN<double> _G1;
    VecN<double> _G2;

    VecN<double> _q;

    double eff_grad(Supercell &cell, VecN<double> &ax, std::unique_ptr<PotentialBase> &ff);
};
