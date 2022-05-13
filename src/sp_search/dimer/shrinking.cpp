#include "sp_search/dimer/shrinking.hpp"

#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

ShrinkingDimer ShrinkingDimer::load(toml::v2::table const &config) {
    ShrinkingDimer opt;

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

bool ShrinkingDimer::find_sp(Supercell &cell, VecN<double> &v, std::unique_ptr<PotentialBase> &ff) {
    double l = _opt.l0;
    double bet = _opt.delta_t;
    double gam = _opt.delta_t;

    _x = cell.activ.view();

    for (size_t i = 0; i < _opt.iter_max_tran; i++) {
        // Compute gradients at end points of dimer
        cell.activ.view() = _x + (0.5 * l) * v;
        ff->gradient(cell, _G1);

        cell.activ.view() = _x - (0.5 * l) * v;
        ff->gradient(cell, _G2);

        using std::swap;
        // Compute d_k (torque)
        swap(_d_p, _d);
        _d = (1 / l) * (_G1 - _G2);
        // double curv = dot(_d, v) / 2;
        _d -= dot(_d, v) * v;

        // Compute g_k (effective force)
        swap(_g_p, _g);
        _g = -0.5 * (_G1 + _G2);
        _g -= 2 * dot(_g, v) * v;

        // Convergence criterion
        double norm_g = norm(_g);

        if (norm_g < _opt.f2norm) {
            cell.activ.view() = _x;
            return dot(_G1 - _G2, v) < 0;
        }

        // Compute beta/gamma
        if (i > 0) {
            switch (_opt.BB_method) {
                case 1:
                    bet = dot(_x - _x_p, _x - _x_p) / dot(_x - _x_p, _g_p - _g);
                    gam = dot(v - _v_p, v - _v_p) / dot(v - _v_p, _d - _d_p);
                    break;
                case 2:
                    bet = dot(_x - _x_p, _g_p - _g) / dot(_g_p - _g, _g_p - _g);
                    gam = dot(v - _v_p, _d - _d_p) / dot(_d - _d_p, _d - _d_p);
                    break;
                default:
                    ALWAYS_CHECK(false, "Invaldid BB method");
            }
        }

        bet = std::min(_opt.s_max / norm_g, std::abs(bet));
        gam = std::abs(gam);

        // dump_supercell(cell, "olkmc.xyz", true);

        // std::cout << std::scientific << std::setprecision(6) << "it: " << i << "\tf " << norm_g
        //           << "\tbet " << bet * norm_g << "\tgam " << gam << "\tl " << l << "\tcurv " <<
        //           curv
        //           << '\n';

        // Update equations of motion
        _x_p = _x;
        _x += bet * _g;

        _v_p = v;
        v -= gam * _d;
        v *= 1 / norm(v);

        l = std::max(_opt.l_min, l / (1 + _opt.delta_t));
    }

    return false;
}
