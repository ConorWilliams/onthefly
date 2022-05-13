#pragma once

#include <bits/c++config.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>

#include "config.hpp"
#include "minimise/LBFGS/core.hpp"
#include "minimise/minimiser_base.hpp"
#include "potentials/potential_base.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

struct DimerRotor {
    std::size_t n = 6;              // LBFGS rotation-history size
    std::size_t iter_max_rot = 20;  // Maximum number of iterations during rotation

    double delta_r;    // (Angstroms) Dimer length
    double theta_tol;  // (Rad) For rotation considered converged

    static DimerRotor load(toml::v2::table const &config);
};

struct Dimer {
    //
    DimerRotor opt;

    std::size_t n = 10;                // Number of previous steps held in memory
    std::size_t iter_max_tran = 2000;  // Number of translations before exit
    std::size_t convex_max = 5;        // Number of consecutive +Ve curvature steps before exit

    bool no_relax_in_convex = false;  // If true, in convex region use only parallel force component

    double boost_parallel = 0;  // Biases effective-force toward parallel component
    double proj_tol = 0;        // Trust tolerance
    double grow_trust = 1.5;    // Trust radius expansion rate
    double shrink_trust = 0.5;  // Trust radius contraction rate

    double f2norm;     // Force convergence criterion (eV/Angstroms)
    double nudge;      // (Angstroms) Distance perturbed along min-mode at SP
    double basin_tol;  // (Angstroms) L2 between active atoms to be considered distinct basins
    double max_trust;  // Maximum trust radius (Angstroms)
    double min_trust;  // Minimum trust radius (Angstroms)

    static Dimer load(toml::v2::table const &config);
};

}  // namespace options

// Class responsible for rotating a dimer (cell+axis) to align axis with minimum-mode
class DimerRotor {
  public:
    DimerRotor(options::DimerRotor const &opt) : _opt(opt), _core(opt.n) {}

    double align(Supercell &cell, VecN<double> &ax, std::unique_ptr<PotentialBase> &ff);

    VecN<double> const &grad() const { return _g0; }

  private:
    options::DimerRotor _opt;

    CoreLBFGS _core;

    VecN<double> _active;  // Store active atoms
    VecN<double> _axisp;   // Temporary axis
    VecN<double> _g0;      // Central grad
    VecN<double> _g1;      // End grad
    VecN<double> _g1p;     // Temp end grad
    VecN<double> _delta_g;
    VecN<double> _theta;
};

// Class responsible for moving a dimer (cell+axis) to a near-by saddle-point
class Dimer {
  public:
    Dimer(options::Dimer const &opt) : _opt(opt), _core(opt.n), _rotor{opt.opt} {}

    bool find_sp(Supercell &cell, VecN<double> &ax, std::unique_ptr<PotentialBase> &ff);

    options::Dimer const &get_opt() const { return _opt; }

  private:
    options::Dimer _opt;

    CoreLBFGS _core;

    DimerRotor _rotor;

    VecN<double> _q;     // Newton step
    VecN<double> _geff;  // Effective gradient

    double eff_grad(Supercell &cell, VecN<double> &ax, std::unique_ptr<PotentialBase> &ff);
};
