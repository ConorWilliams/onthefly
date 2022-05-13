#include "sp_search/dimer/dimer.hpp"

#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

DimerRotor DimerRotor::load(toml::v2::table const &config) {
    DimerRotor opt;

    opt.n = config["sp_search"]["dimer"]["n_rot"].value_or(opt.n);
    opt.iter_max_rot = config["sp_search"]["dimer"]["iter_max_rot"].value_or(opt.iter_max_rot);

    opt.delta_r = fetch<double>(config, "sp_search", "dimer", "delta_r");
    opt.theta_tol = fetch<double>(config, "sp_search", "dimer", "theta_tol");

    return opt;
}

Dimer Dimer::load(toml::v2::table const &config) {
    Dimer opt;

    opt.opt = DimerRotor::load(config);

    opt.n = config["sp_search"]["dimer"]["n_tran"].value_or(opt.n);
    opt.iter_max_tran = config["sp_search"]["dimer"]["iter_max_tran"].value_or(opt.iter_max_tran);
    opt.convex_max = config["sp_search"]["dimer"]["convex_max"].value_or(opt.convex_max);

    opt.no_relax_in_convex
        = config["sp_search"]["no_relax_in_convex"]["no_relax_in_convex"].value_or(
            opt.no_relax_in_convex);

    opt.boost_parallel
        = config["sp_search"]["dimer"]["boost_parallel"].value_or(opt.boost_parallel);
    opt.proj_tol = config["sp_search"]["dimer"]["proj_tol"].value_or(opt.proj_tol);
    opt.grow_trust = config["sp_search"]["dimer"]["grow_trust"].value_or(opt.grow_trust);
    opt.shrink_trust = config["sp_search"]["dimer"]["shrink_trust"].value_or(opt.shrink_trust);

    opt.f2norm = fetch<double>(config, "sp_search", "dimer", "f2norm");
    opt.nudge = fetch<double>(config, "sp_search", "dimer", "nudge");
    opt.basin_tol = fetch<double>(config, "sp_search", "dimer", "basin_tol");
    opt.max_trust = fetch<double>(config, "sp_search", "dimer", "max_trust");
    opt.min_trust = fetch<double>(config, "sp_search", "dimer", "min_trust");

    return opt;
}

}  // namespace options

// Rotate the dimer axis to align with the minimum eigen-mode of the Hessian, cell is modified but
// returned as recived. After calling this function get_grad() will return the gradient at the
// centre of the dimer
double DimerRotor::align(Supercell &cell, VecN<double> &ax, std::unique_ptr<PotentialBase> &ff) {
    //
    _active = cell.activ.view();  // Save active
    _core.clear();

    ff->gradient(cell, _g0);  // Gradient at centre

    cell.activ.view() = _active + _opt.delta_r * ax;
    ff->gradient(cell, _g1);

    for (size_t i = 0;; i++) {
        _delta_g = _g1 - _g0;
        _delta_g -= dot(_delta_g, ax) * ax;  // Torque

        _core(ax, _delta_g, _theta);  // Use lbfgs to find rotation plane

        _theta -= dot(_theta, ax) * ax;  // Ortho
        _theta *= 1 / norm(_theta);      //      normalization

        double b_1 = dot(_g1 - _g0, _theta) / _opt.delta_r;
        double c_x0 = dot(_g1 - _g0, ax) / _opt.delta_r;
        double theta_1 = -0.5 * std::atan(b_1 / std::abs(c_x0));  // Trial rotation angle

        if (std::abs(theta_1) < _opt.theta_tol || i == _opt.iter_max_rot) {
            cell.activ.view() = _active;  // Return cell to original state
            return c_x0;
        } else {
            _axisp = ax * std::cos(theta_1) + _theta * std::sin(theta_1);

            cell.activ.view() = _active + _opt.delta_r * _axisp;
            ff->gradient(cell, _g1p);

            double c_x1 = dot(_g1p - _g0, _axisp) / _opt.delta_r;
            double a_1 = (c_x0 - c_x1 + b_1 * sin(2 * theta_1)) / (1 - std::cos(2 * theta_1));
            double theta_min = 0.5 * std::atan(b_1 / a_1);  // Optimal rotation

            // Flip if extrema is maxima
            if (a_1 * std::cos(2 * theta_min) - a_1 + b_1 * std::sin(2 * theta_min) > 0) {
                theta_min += M_PI / 2;
            }

            // std::cout << i << " theta " << theta_min << " curv " << c_x0 << std::endl;

            ax = ax * std::cos(theta_min) + _theta * std::sin(theta_min);

            // Interpolate force at new rotation
            _g1 = (std::sin(theta_1 - theta_min) / std::sin(theta_1)) * _g1
                  + (std::sin(theta_min) / std::sin(theta_1)) * _g1p
                  + (1 - std::cos(theta_min) - std::sin(theta_min) * std::tan(0.5 * theta_1)) * _g0;

            if (std::abs(theta_min) < _opt.theta_tol) {
                cell.activ.view() = _active;  // Return cell to original state
                return dot(_g1 - _g0, ax) / _opt.delta_r;
            }
        }
    }
}

bool Dimer::find_sp(Supercell &cell, VecN<double> &ax, std::unique_ptr<PotentialBase> &ff) {
    // Uses LBFGS loop with modified curvature escape clauses

    _core.clear();

    double curv = eff_grad(cell, ax, ff);

    double trust = _opt.min_trust;

    std::size_t convex_count = 0;

    for (std::size_t i = 0; i < _opt.iter_max_tran; ++i) {
        // dump_supercell(cell, "olkmc.xyz", true);

        // std::cout << std::scientific << std::setprecision(6) << "it: " << i << "\tf "
        //           << dot(_rotor.grad(), _rotor.grad()) << "\ttr " << trust << "\tcv " << curv
        //           << '\n';

        if (dot(_geff, _geff) < _opt.f2norm * _opt.f2norm) {
            return true;
        } else if (convex_count >= _opt.convex_max) {
            return false;
        }

        _core(cell.activ.view(), _geff, _q);

        // Trust-radius based line-search;
        cell.activ.view() -= std::min(1.0, trust / norm(_q)) * _q;

        curv = eff_grad(cell, ax, ff);

        if (double proj = dot(_geff, _q); proj < -_opt.proj_tol) {
            trust = std::max(_opt.min_trust, _opt.shrink_trust * trust);
        } else if (proj > _opt.proj_tol) {
            trust = std::min(_opt.max_trust, _opt.grow_trust * trust);
        }

        if (curv > 0) {
            convex_count += 1;
        } else {
            convex_count = 0;
        }
    }

    // std::cerr << "Failed to converge dimer!" << std::endl;

    return false;
}

// Transform the gradient at the centre into the modified gradient
double Dimer::eff_grad(Supercell &cell, VecN<double> &ax, std::unique_ptr<PotentialBase> &ff) {
    //

    double curv = _rotor.align(cell, ax, ff);  // cell.activ is overwritten

    double mag = norm(_rotor.grad());

    if (_opt.no_relax_in_convex && curv > 0) {
        _geff = -dot(_rotor.grad(), ax) * ax;
    } else {
        _geff = _rotor.grad() - (2 + _opt.boost_parallel) * dot(_rotor.grad(), ax) * ax;
    }

    _geff *= mag / norm(_geff);

    return curv;
}
