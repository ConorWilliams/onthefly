#pragma once

#include "config.hpp"
#include "minimise/LBFGS/core.hpp"
#include "minimise/minimiser_base.hpp"
#include "supercell.hpp"
#include "utility.hpp"

namespace options {

struct MinimiseHYBRID {
    std::size_t iter_max_tran = 1000;  // Number of translations before exit
    std::size_t BB_method = 1;         // Use BB1 vs BB2 https://doi.org/10.1137/19M1253356

    double delta_t;  // Shrinking rate
    double l0;       // (Angstroms) Dimer length
    double l_min;    // (Angstrom) Min dimer length
    double f2norm;   // Force convergence criterion (eV/Angstroms)
    double boost;    // Multiplier by which gradient is boosted in negative curvature regions

    std::size_t n = 10;         // Number of previous steps held in memory
    double max_trust = 0.5;     // Maximum trust radius (Angstroms)
    double min_trust = 0.05;    // Minimum trust radius (Angstroms)
    double grow_trust = 1.5;    // Trust radius expansion rate
    double shrink_trust = 0.5;  // Trust radius contraction rate

    static MinimiseHYBRID load(toml::v2::table const &config);
};

}  // namespace options

// Minimiser that uses the LBFGS algorithm for dimer translation and OSD-BB method for dimer
// rotation
class MinimiseHYBRID final : public MinimiserBase {
  public:
    std::unique_ptr<MinimiserBase> clone() const override;

    MinimiseHYBRID(options::MinimiseHYBRID const &opt) : _opt(opt), _core(_opt.n) {}

    bool minimise(Supercell &, std::unique_ptr<PotentialBase> &, VecN<double> const &) override;

    bool minimise(Supercell &cell, std::unique_ptr<PotentialBase> &ff) override {
        VecN<double> rand(cell.activ.size() * 3);
        rand_normal(rand);
        return minimise(cell, ff, rand);
    }

  private:
    options::MinimiseHYBRID _opt;

    CoreLBFGS _core;

    VecN<double> _x;
    VecN<double> _v;

    VecN<double> _G1;
    VecN<double> _G2;

    VecN<double> _g;
    VecN<double> _d;

    VecN<double> _d_p;
    VecN<double> _v_p;

    VecN<double> _q;
};
