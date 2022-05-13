#include "kinetics/basin.hpp"

#include <random>
#include <string>

#include "local/geometry.hpp"
#include "utility.hpp"

namespace {

constexpr double invBoltz = 16021766340.0 / 1380649.0;  // eV^{-1}

}

namespace options {

Basin Basin::load(toml::v2::table const &config) {
    Basin opt;

    opt.temp = fetch<double>(config, "kinetics", "temperature");
    opt.max_barrier = config["kinetics"]["max_barrier"].value_or(opt.max_barrier);

    return opt;
}

}  // namespace options

Mechanism const &Basin::local_mech::onto(Supercell &cell, std::vector<Geometry> &geos) const {
    //
    CHECK(_atom_idx < cell.activ.size(), "Invalid atom index");

    // std::cout << "abs " << env[_atom_idx]->mechs[_mech_off].abs_cap << std::endl;

    // permute required if cell changed due to switcheroo(superbasin)
    std::optional<Geometry::Result> Res = geos[_atom_idx].permute_onto(_env->delta, _env->geo);

    ALWAYS_CHECK(Res, "unable to align cell & mechs geos");

    Res->R.transposeInPlace();

    std::size_t j = 0;

    for (std::size_t i = 0; i < geos[_atom_idx].size(); i++) {
        if (geos[_atom_idx][i].col.state == Colour::activ) {
            cell.activ[geos[_atom_idx][i].idx].vec
                += (Res->R * _env->mechs[_mech_off].disp[j++].matrix()).array();
        }
    }

    return _env->mechs[_mech_off];
}

void Basin::local_mech::refine(std::vector<Geometry> &geo) const { _env.refine(geo[_atom_idx]); }

Basin::Basin(options::Basin const &opt,
             Supercell const &cell,
             std::vector<Catalogue::pointer> const &env)
    : _state(cell.activ.view()) {
    for (std::size_t i = 0; i < env.size(); ++i) {
        //
        // std::cout << env[i]->mechs.size() << std::endl;
        // if (env[i]->mechs.size() > 100) {
        //     auto m0 = env[i]->mechs[0];

        //     for (auto const &m : env[i]->mechs) {
        //         double sum_sq = 0;

        //         for (std::size_t i = 0; i < m.disp.size(); i++) {
        //             sum_sq += norm_sq(m0.disp[i] - m.disp[i]);
        //         }

        //         std::cout << "sum: " << std::sqrt(sum_sq) << std::endl;

        //     }

        //     exit(1);
        // }

        for (size_t j = 0; j < env[i]->mechs.size(); j++) {
            double fwd = env[i]->mechs[j].activ_energy;

            if (fwd < opt.max_barrier) {
                CHECK(fwd > 0, "Negative energy barrier!");

                double rate = env[i]->mechs[j].pre_factor * std::exp(fwd / opt.temp * -invBoltz);

                double rev = fwd - env[i]->mechs[j].delta_energy;

                // if (fwd < 1e-3 || rev < 1e-3) {
                //     // TODO: remove hack, this prevents sim getting stuck in ledge mechs
                //     continue;
                // }

                // if (true) {
                //     std::cout << fwd << std::endl;
                // }

                _mechs.emplace_back(rate, std::max(fwd, rev), i, j, env[i]);

                _rate_sum += rate;
            }
        }
    }
}

Basin::choice Basin::kmc_choice(pcg64 &psudo_rng, std::size_t basin) const {
    std::uniform_real_distribution<double> uniform;

    std::size_t const mech = [&]() {
        double const lim = uniform(psudo_rng) * rate_sum();
        double sum = 0;

        for (std::size_t i = 0; i < size(); ++i) {
            if (_mechs[i].exit_mech) {
                sum += _mechs[i].rate;
                if (sum > lim) {
                    return i;
                }
            }
        }

        throw std::runtime_error("Hit end of normal choice: " + std::to_string(size()));
    }();

    double const rate = _mechs[mech].rate;

    std::cout << "Rate: " << rate << " : " << rate / rate_sum() * 100 << "% of ";
    std::cout << size() << " choices\n";

    return {
        false,
        mech,
        -std::log(uniform(psudo_rng)) / _rate_sum,
        basin,
    };
}
