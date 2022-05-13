#include "sp_search/find_mech.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <stdexcept>

#include "config.hpp"
#include "local/environment.hpp"
#include "pcg_random.hpp"
#include "sp_search/vineyard.hpp"
#include "supercell.hpp"
#include "utility.hpp"

namespace options {

FindMechanisms FindMechanisms::load(toml::v2::table const& config) {
    FindMechanisms opt;

    opt.consecutive = config["sp_search"]["consecutive"].value_or(opt.consecutive);
    opt.max_search = config["sp_search"]["max_search"].value_or(opt.max_search);
    opt.vineyard = config["sp_search"]["vineyard"].value_or(opt.vineyard);
    opt.vine_zero_tol = config["sp_search"]["vine_zero_tol"].value_or(opt.vine_zero_tol);

    opt.r_perturbation = fetch<double>(config, "sp_search", "r_perturbation");
    opt.stddev = fetch<double>(config, "sp_search", "stddev");

    opt.proto = Mechanism::load(config);

    if (!opt.vineyard) {
        opt.const_pre_factor = fetch<double>(config, "sp_search", "const_pre_factor");
    }

    return opt;
}

}  // namespace options

namespace {

void random_local_pertubation(std::size_t n, Supercell& dimer, double range, double stddev) {
    // Seed with a real random value, if available
    static thread_local pcg64 rng(pcg_extras::seed_seq_from<std::random_device>{});

    static thread_local std::normal_distribution<> mag_dist(1, 1);

    static thread_local std::normal_distribution<> normal_dist;

    static thread_local std::normal_distribution<> guassian(0, stddev);

    Vec3<double> centre = dimer.activ[n].vec;

    double mag = mag_dist(rng);

    for (std::size_t i = 0; i < dimer.activ.size() * 3; i += 3) {
        //
        Vec3<double> delta = dimer.min_image(dimer.activ[i / 3].vec, centre);

        if (norm_sq(delta) > range * range) {
            continue;
        }

        Vec3<double> dr = {
            normal_dist(rng),
            normal_dist(rng),
            normal_dist(rng),
        };

        dr /= norm(dr);

        // dr is now unit vec with random direction

        double cut = std::exp(-norm_sq(delta) / (range * range));

        dimer.activ.view().segment<3>(i) += mag * cut * guassian(rng) * dr;
    }

    return;
}

// Compute change in active atoms positions, attempts to correct for com drift
VecN<double> mech_disp(Supercell const& xi, Supercell const& xf) {
    CHECK(xi.activ.size() == xf.activ.size(), "Number of atoms are different");

    VecN<double> dr = xf.activ.view() - xi.activ.view();

    // Can skip translational correction if bound atoms are present as they prevent drift
    if (xi.bound.size() > 0) {
        return dr;
    }

    // Find the atom with max atomic displacement
    std::size_t j = 0;

    {
        double max = norm_sq(dr.segment<3>(0));

        for (int i = 1; i < dr.size() / 3; ++i) {
            double ns = norm_sq(dr.segment<3>(3 * i));
            if (ns > max) {
                max = ns;
                j = i;
            }
        }
    }

    // Compute the com of all atoms far from the central atom
    std::size_t count = 0;
    Vec3<double> com_i = Vec3<double>::Zero();
    Vec3<double> com_f = Vec3<double>::Zero();

    for (size_t i = 0; i < xi.activ.size(); i++) {
        double ni = norm_sq(xi.min_image(xi.activ[i].vec, xi.activ[j].vec));
        double nf = norm_sq(xf.min_image(xf.activ[i].vec, xf.activ[j].vec));

        if (ni > 6 * 6 && nf > 6 * 6) {
            com_i += xi.activ[i].vec;
            com_f += xf.activ[i].vec;
            ++count;
        }
    }

    ALWAYS_CHECK(count > 0, "no atoms :(");

    // Compute com drift
    VecN<double> d_com = (com_f - com_i) / count;

    // Apply com-drift correction
    for (int i = 0; i < dr.size() / 3; ++i) {
        dr.segment<3>(3 * i) -= d_com;
    }

    return dr;
}

}  // namespace

std::vector<ProtoMech> find_mechanisms(options::FindMechanisms const& opt,
                                       Workcell const& init,
                                       std::unique_ptr<PotentialBase>& ff,
                                       std::unique_ptr<SearchBase>& finder) {
    // New mechs stored here
    std::vector<ProtoMech> mechs;
    //
    Supercell dimer = init;
    Supercell final = init;

    std::size_t count = 0;  // count consecutive failure/mech rediscoveries

    Vineyard vine(opt.vine_zero_tol);  // for computing harmonic prefactor

    // Vineyard vine2(opt.vine_zero_tol);

    if (opt.vineyard) {
        vine.load_basin(init, ff);
    }

    for (std::size_t i = 0; count < opt.consecutive; i++) {
        if (i >= opt.max_search) {
            // std::cout << "WARNING: find_mechanisms hit max_search, consider increasing it\n";
            break;
        }

        dimer = init;

        random_local_pertubation(init.centre, dimer, opt.r_perturbation, opt.stddev);

        count++;

        try {
            if (finder->find_sp(init, dimer, final, ff)) {
                double Ei = ff->energy(init);
                double Es = ff->energy(dimer);
                double Ef = ff->energy(final);

                ProtoMech mech{{Es - Ei, Ef - Ei, opt.const_pre_factor}, mech_disp(init, final)};

                auto it = std::find_if(mechs.begin(), mechs.end(), [&](ProtoMech const& other) {
                    return mech.within_tol(other,
                                           opt.proto.energy_abs_tol,
                                           opt.proto.energy_frac_tol,
                                           opt.proto.r_tol);
                });

                if (it == mechs.end()) {
                    if (opt.vineyard) {
                        if (vine.load_sp(dimer, ff)) {
                            mech.pre_factor = vine.pre_factor();
                            mechs.push_back(std::move(mech));
                            count = 0;
                        }
                    } else {
                        mechs.push_back(std::move(mech));
                        count = 0;
                    }
                }
            }
        } catch (std::runtime_error const& err) {
            std::cerr << "Caught: " << err.what() << std::endl;
            count--;
        }
    }

    return mechs;
}
