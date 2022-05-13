#pragma once

#include "config.hpp"
#include "minimise/minimiser_base.hpp"
#include "supercell.hpp"
#include "utility.hpp"

namespace options {

struct MinimiseBB {
    std::size_t iter_max_tran = 1000;  // Number of translations before exit
    std::size_t BB_method = 1;         // Use BB1 vs BB2 https://doi.org/10.1137/19M1253356

    double delta_t;  // Shrinking rate
    double l0;       // (Angstroms) Dimer length
    double l_min;    // (Angstrom) Min dimer length
    double s_max;    // (Angstrom) max step size
    double f2norm;   // Force convergence criterion (eV/Angstroms)

    double nudge;      // (Angstroms) Distance perturbed along min-mode at SP
    double basin_tol;  // (Angstroms) L2 between active atoms to be considered distinct basins

    static MinimiseBB load(toml::v2::table const &config);
};

}  // namespace options

// Minimiser that uses the OSD-BB algorithm for dimer translation and dimer
// rotation
class MinimiseBB final : public MinimiserBase {
  public:
    std::unique_ptr<MinimiserBase> clone() const override;

    MinimiseBB(options::MinimiseBB const &opt) : _opt(opt) {}

    bool minimise(Supercell &, std::unique_ptr<PotentialBase> &, VecN<double> const &) override;

    bool minimise(Supercell &cell, std::unique_ptr<PotentialBase> &ff) override {
        VecN<double> rand(cell.activ.size() * 3);
        rand_normal(rand);
        return minimise(cell, ff, rand);
    }

  private:
    options::MinimiseBB _opt;

    VecN<double> _x;
    VecN<double> _v;

    VecN<double> _G1;
    VecN<double> _G2;

    VecN<double> _g;
    VecN<double> _d;

    VecN<double> _x_p;
    VecN<double> _g_p;
    VecN<double> _d_p;
    VecN<double> _v_p;
};
