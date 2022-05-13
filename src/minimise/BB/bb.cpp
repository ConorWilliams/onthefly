#include "minimise/BB/bb.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

MinimiseBB MinimiseBB::load(toml::v2::table const &config) {
    MinimiseBB opt;

    opt.iter_max_tran = config["minimiser"]["mod"]["iter_max_tran"].value_or(opt.iter_max_tran);
    opt.BB_method = config["minimiser"]["mod"]["BB_method"].value_or(opt.BB_method);

    opt.l0 = fetch<double>(config, "minimiser", "mod", "l0");
    opt.l_min = fetch<double>(config, "minimiser", "mod", "l_min");
    opt.s_max = fetch<double>(config, "minimiser", "mod", "s_max");
    opt.delta_t = fetch<double>(config, "minimiser", "mod", "delta_t");
    opt.f2norm = fetch<double>(config, "minimiser", "mod", "f2norm");
    opt.nudge = fetch<double>(config, "minimiser", "mod", "nudge");
    opt.basin_tol = fetch<double>(config, "minimiser", "mod", "basin_tol");

    return opt;
}

}  // namespace options

std::unique_ptr<MinimiserBase> MinimiseBB::clone() const {
    return std::make_unique<MinimiseBB>(*this);
}

bool MinimiseBB::minimise(Supercell &cell,
                          std::unique_ptr<PotentialBase> &ff,
                          VecN<double> const &hint_v) {
    //
    double l = _opt.l0;
    double bet = _opt.delta_t;
    double gam = _opt.delta_t;

    _x = cell.activ.view();
    _v = hint_v;

    for (size_t i = 0; i < _opt.iter_max_tran; i++) {
        // Compute gradients at end points of dimer
        cell.activ.view() = _x + (0.5 * l) * _v;
        ff->gradient(cell, _G1);

        cell.activ.view() = _x - (0.5 * l) * _v;
        ff->gradient(cell, _G2);

        using std::swap;

        // Compute g_k (effective force)
        swap(_g_p, _g);
        _g = -0.5 * (_G1 + _G2);

        // Convergence criterion
        double norm_g = norm(_g);

        // Compute d_k
        swap(_d_p, _d);
        _d = (1 / l) * (_G1 - _G2);
        double curv = 0.5 * dot(_d, _v);
        _d -= dot(_d, _v) * _v;

        // Boost gradient if in flat neg curve region
        if (curv < 0) {
            _g += (_opt.f2norm * 100 / norm_g) * _g;
            norm_g = norm(_g);
        }

        if (norm_g < _opt.f2norm) {
            cell.activ.view() = _x;
            return true;
        }

        // Compute beta/gamma
        if (i > 0) {
            switch (_opt.BB_method) {
                case 1:
                    bet = dot(_x - _x_p, _x - _x_p) / dot(_x - _x_p, _g_p - _g);
                    gam = dot(_v - _v_p, _v - _v_p) / dot(_v - _v_p, _d - _d_p);
                    break;
                case 2:
                    bet = dot(_x - _x_p, _g_p - _g) / dot(_g_p - _g, _g_p - _g);
                    gam = dot(_v - _v_p, _d - _d_p) / dot(_d - _d_p, _d - _d_p);
                    break;
                default:
                    ALWAYS_CHECK(false, "Invaldid BB method");
            }
        }

        bet = std::min(_opt.s_max / norm_g, std::abs(bet));
        gam = std::abs(gam);

        dump_supercell(cell, "olkmc.xyz", true);

        std::cout << std::scientific << std::setprecision(6) << "it: " << i << "\tf " << norm_g
                  << "\tbet " << bet * norm_g << "\tgam " << gam << "\tl " << l << "\tcurv " << curv
                  << '\n';

        // Update equations of motion
        _x_p = _x;
        _x += bet * _g;

        _v_p = _v;
        _v -= gam * _d;
        _v *= 1 / norm(_v);

        l = std::max(_opt.l_min, l / (1 + _opt.delta_t));
    }

    return false;
}