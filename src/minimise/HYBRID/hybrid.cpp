#include "minimise/HYBRID/hybrid.hpp"

#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

MinimiseHYBRID MinimiseHYBRID::load(toml::v2::table const &config) {
    MinimiseHYBRID opt;

    opt.iter_max_tran = config["minimiser"]["hybrid"]["iter_max_tran"].value_or(opt.iter_max_tran);
    opt.BB_method = config["minimiser"]["hybrid"]["BB_method"].value_or(opt.BB_method);

    opt.l0 = fetch<double>(config, "minimiser", "hybrid", "l0");
    opt.l_min = fetch<double>(config, "minimiser", "hybrid", "l_min");
    opt.delta_t = fetch<double>(config, "minimiser", "hybrid", "delta_t");
    opt.f2norm = fetch<double>(config, "minimiser", "hybrid", "f2norm");
    opt.boost = fetch<double>(config, "minimiser", "hybrid", "boost");

    opt.max_trust = config["minimiser"]["hybrid"]["max_trust"].value_or(opt.max_trust);
    opt.min_trust = config["minimiser"]["hybrid"]["min_trust"].value_or(opt.min_trust);
    opt.grow_trust = config["minimiser"]["hybrid"]["grow_trust"].value_or(opt.grow_trust);
    opt.shrink_trust = config["minimiser"]["hybrid"]["shrink_trust"].value_or(opt.shrink_trust);
    opt.n = config["minimiser"]["hybrid"]["n"].value_or(opt.n);

    return opt;
}

}  // namespace options

std::unique_ptr<MinimiserBase> MinimiseHYBRID::clone() const {
    return std::make_unique<MinimiseHYBRID>(*this);
}

bool MinimiseHYBRID::minimise(Supercell &cell,
                              std::unique_ptr<PotentialBase> &ff,
                              VecN<double> const &hint_v) {
    double l = _opt.l0;
    double gam = _opt.delta_t;
    double trust = _opt.min_trust;

    _x = cell.activ.view();
    _v = hint_v;

    for (size_t i = 0; i < _opt.iter_max_tran; i++) {
        // Compute gradients at end points of dimer
        cell.activ.view() = _x + (0.5 * l) * _v;
        ff->gradient(cell, _G1);

        cell.activ.view() = _x - (0.5 * l) * _v;
        ff->gradient(cell, _G2);

        // Compute g_k (gradient at centre)
        _g = 0.5 * (_G1 + _G2);

        // Adjust trust radius (max step size)
        if (i > 0) {
            if (dot(_g, _q) < 0) {
                trust = std::max(_opt.min_trust, _opt.shrink_trust * trust);
            } else {
                trust = std::min(_opt.max_trust, _opt.grow_trust * trust);
            }
        }

        // Compute d_k
        using std::swap;
        swap(_d_p, _d);
        _d = (1 / l) * (_G1 - _G2);
        double curv = 0.5 * dot(_d, _v);
        _d -= dot(_d, _v) * _v;

        // Boost gradient if in neg curve region
        if (curv < 0) {
            _g += (_opt.f2norm * _opt.boost / norm(_g)) * _g;
        }

        double norm_g = norm(_g);

        // Convergence criterion == strict minima
        if (norm_g < _opt.f2norm && curv > 0) {
            cell.activ.view() = _x;
            return true;
        }

        // Compute beta/gamma
        if (i > 0) {
            switch (_opt.BB_method) {
                case 1:
                    gam = dot(_v - _v_p, _v - _v_p) / dot(_v - _v_p, _d - _d_p);
                    break;
                case 2:
                    gam = dot(_v - _v_p, _d - _d_p) / dot(_d - _d_p, _d - _d_p);
                    break;
                default:
                    ALWAYS_CHECK(false, "Invaldid BB method");
            }

            gam = std::abs(gam);
        }

        // Get descent directrion via lbfgs
        _core(_x, _g, _q);

        // dump_supercell(cell, "olkmc.xyz", true);

        // std::cout << std::scientific << std::setprecision(6) << "it: " << i << "\tf " << norm_g
        //           << "\tgam " << gam << "\tl " << l << "\tcurv " << curv << "\ttrust " << trust
        //           << '\n';

        // Update equations of motion
        _x -= std::min(1.0, trust / norm(_q)) * _q;

        _v_p = _v;
        _v -= gam * _d;
        _v *= 1 / norm(_v);

        l = std::max(_opt.l_min, l / (1 + _opt.delta_t));

        ////
    }

    return false;
}