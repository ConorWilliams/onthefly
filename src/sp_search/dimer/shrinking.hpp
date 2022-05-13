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

struct ShrinkingDimer {
    std::size_t iter_max_tran = 1000;  // Number of translations before exit
    std::size_t BB_method = 1;         // Use BB1 vs BB2 https://doi.org/10.1137/19M1253356

    double delta_t;  // Shrinking rate
    double l0;       // (Angstroms) Dimer length
    double l_min;    // (Angstrom) Min dimer length
    double s_max;    // (Angstrom) max step size
    double f2norm;   // Force convergence criterion (eV/Angstroms)

    double nudge;      // (Angstroms) Distance perturbed along min-mode at SP
    double basin_tol;  // (Angstroms) L2 between active atoms to be considered distinct basins

    static ShrinkingDimer load(toml::v2::table const &config);
};

}  // namespace options

// Class responsible for moving a dimer (cell+axis) to a near-by saddle-point
class ShrinkingDimer {
  public:
    ShrinkingDimer(options::ShrinkingDimer const &opt) : _opt(opt) {}

    bool find_sp(Supercell &cell, VecN<double> &ax, std::unique_ptr<PotentialBase> &ff);

    options::ShrinkingDimer const &get_opt() const { return _opt; }

  private:
    options::ShrinkingDimer _opt;

    VecN<double> _x;

    VecN<double> _G1;
    VecN<double> _G2;

    VecN<double> _g;
    VecN<double> _d;

    VecN<double> _x_p;
    VecN<double> _g_p;
    VecN<double> _d_p;
    VecN<double> _v_p;

    double eff_grad(Supercell &cell, VecN<double> &ax, std::unique_ptr<PotentialBase> &ff);
};
