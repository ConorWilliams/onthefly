#include "sp_search/dimer/l_shrink.hpp"

#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

LShrink LShrink::load(toml::v2::table const &config) {
    LShrink opt;

    opt.iter_max_tran = config["sp_search"]["shrink"]["iter_max_tran"].value_or(opt.iter_max_tran);
    opt.BB_method = config["sp_search"]["shrink"]["BB_method"].value_or(opt.BB_method);

    opt.l0 = fetch<double>(config, "sp_search", "shrink", "l0");
    opt.l_min = fetch<double>(config, "sp_search", "shrink", "l_min");
    opt.s_max = fetch<double>(config, "sp_search", "shrink", "s_max");
    opt.delta_t = fetch<double>(config, "sp_search", "shrink", "delta_t");
    opt.f2norm = fetch<double>(config, "sp_search", "shrink", "f2norm");
    opt.nudge = fetch<double>(config, "sp_search", "shrink", "nudge");
    opt.basin_tol = fetch<double>(config, "sp_search", "shrink", "basin_tol");

    return opt;
}

}  // namespace options

bool LShrink::find_sp(Supercell &cell, VecN<double> &v, std::unique_ptr<PotentialBase> &ff) {
    //
    auto const n = v.size();
    _x.resize(2 * n);
    _gx.resize(2 * n);

    _x.head(n) = cell.activ.view();
    _x.tail(n) = v;

    for (size_t i = 0; i < _opt.iter_max_tran; i++) {
        // Compute gradients at end points of dimer
        cell.activ.view() = _x.head(n) + (0.5 * _opt.l0) * _x.tail(n);
        ff->gradient(cell, _G1);

        cell.activ.view() = _x.head(n) - (0.5 * _opt.l0) * _x.tail(n);
        ff->gradient(cell, _G2);

        // Compute d_k
        _gx.tail(n) = (1 / _opt.l0) * (_G1 - _G2);
        _gx.tail(n) -= dot(_gx.tail(n), _x.tail(n)) * _x.tail(n);  // Perpendiularise

        // Compute g_k (effective gradient)
        _gx.head(n) = 0.5 * (_G1 + _G2);
        _gx.head(n) -= 2 * dot(_gx.head(n), _x.tail(n)) * _x.tail(n);

        // Convergence criterion
        double norm_g = norm(_gx.head(n));

        if (norm_g < _opt.f2norm) {
            exit(1);
            cell.activ.view() = _x.head(n);
            v = _x.tail(n);
            return true;
        }

        // Update equations of motion
        _core(_x, _gx, _q);

        //
        CHECK(dot(_gx, _q) > 0, "Ascent direction");

        // dump_supercell(cell, "olkmc.xyz", true);
        // double curv = dot(_G1 - _G2, v) / (2 * _opt.l0);
        // std::cout << std::scientific << std::setprecision(6) << "it: " << i << "\tf " << norm_g
        //           << "\tcurv " << curv << "\tstep " << norm(_q) << '\n';

        _x -= std::min(1.0, 0.1 / norm(_q)) * _q;

        // Renormailisation
        _x.tail(n) *= 1 / norm(_x.tail(n));
    }

    return false;
}
