#include "minimise/LBFGS/lbfgs.hpp"

#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

MinimiseLBFGS MinimiseLBFGS::load(toml::v2::table const &config) {
    MinimiseLBFGS opt;

    opt.n = config["minimiser"]["lbfgs"]["n"].value_or(opt.n);
    opt.iter_max = config["minimiser"]["lbfgs"]["iter_max"].value_or(opt.iter_max);
    opt.f2norm = config["minimiser"]["lbfgs"]["f2norm"].value_or(opt.f2norm);
    opt.proj_tol = config["minimiser"]["lbfgs"]["proj_tol"].value_or(opt.proj_tol);
    opt.max_trust = config["minimiser"]["lbfgs"]["max_trust"].value_or(opt.max_trust);
    opt.min_trust = config["minimiser"]["lbfgs"]["min_trust"].value_or(opt.min_trust);
    opt.grow_trust = config["minimiser"]["lbfgs"]["grow_trust"].value_or(opt.grow_trust);
    opt.shrink_trust = config["minimiser"]["lbfgs"]["shrink_trust"].value_or(opt.shrink_trust);

    return opt;
}

}  // namespace options

std::unique_ptr<MinimiserBase> MinimiseLBFGS::clone() const {
    return std::make_unique<MinimiseLBFGS>(*this);
}

bool MinimiseLBFGS::minimise(Supercell &cell, std::unique_ptr<PotentialBase> &ff) {
    _core.clear();

    ff->gradient(cell, _gx);

    double trust = _opt.min_trust;

    for (std::size_t i = 0; i < _opt.iter_max; ++i) {
        // dump_supercell(cell, "olkmc.xyz", true);

        // std::cout << std::scientific << std::setprecision(6) << "it: " << i << "\tf "
        //           << dot(_gx, _gx) << "\ttr " << trust << '\n';

        if (dot(_gx, _gx) < _opt.f2norm * _opt.f2norm) {
            return true;
        }

        _core(cell.activ.view(), _gx, _q);

        //
        CHECK(dot(_gx, _q) > 0, "Ascent direction");

        // Trust-radius based line-search;
        cell.activ.view() -= std::min(1.0, trust / norm(_q)) * _q;

        ff->gradient(cell, _gx);

        double proj = dot(_gx, _q);

        if (proj < -_opt.proj_tol) {
            trust = std::max(_opt.min_trust, _opt.shrink_trust * trust);
        } else if (proj > _opt.proj_tol) {
            trust = std::min(_opt.max_trust, _opt.grow_trust * trust);
        }
    }

    // std::cerr << "Failed to converge!" << std::endl;

    return false;
}
