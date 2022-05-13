#include <exception>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "local/catalogue.hpp"
#include "local/classify.hpp"
#include "local/environment.hpp"
#include "minimise/minimiser_base.hpp"
#include "olkmc.hpp"
#include "package/package.hpp"
#include "potentials/EAM/data.hpp"
#include "potentials/EAM/potential.hpp"
#include "potentials/adapt.hpp"
#include "potentials/potential_base.hpp"
#include "riften/thiefpool.hpp"
#include "sp_search/dimer/dimer.hpp"
#include "sp_search/find_mech.hpp"
#include "sp_search/sp_search_base.hpp"
#include "sp_search/vineyard.hpp"
#include "stream.hpp"
#include "structopt/app.hpp"
#include "supercell.hpp"
#include "utility.hpp"

struct CommandLineArgs {
    std::string config_file;

    static CommandLineArgs parse(int argc, char *argv[]) try {
        return structopt::app("olkmc").parse<CommandLineArgs>(argc, argv);
    } catch (structopt::exception &e) {
        std::cout << e.what() << "\n";
        std::cout << e.help();
        throw e;
    }
};

STRUCTOPT(CommandLineArgs, config_file);

void update_catalogue(Classify &classify,
                      Catalogue &cat,
                      Packager &packager,
                      Supercell const &cell,
                      std::unique_ptr<SearchBase> const &finder,
                      std::unique_ptr<PotentialBase> const &ff,
                      options::FindMechanisms const &opt_find,
                      std::vector<DiscreteKey> &keys,
                      std::vector<Geometry> &geos,
                      std::vector<Catalogue::pointer> &env) {
    //
    static riften::Thiefpool pool;

    classify(cell, keys, geos);

    if (std::vector cens = cat.canon_update(keys, geos, env); !cens.empty()) {
        //
        std::vector pkgs = packager.pack(cell, cens);

        Bar bar(pkgs.size());

        for (auto &&pk : pkgs) {
            pk.f_mechs = pool.enqueue([&, f = finder->clone(), p = ff->clone()]() mutable {
                auto ret = find_mechanisms(opt_find, pk.subcell, p, f);
                bar.tick();
                return ret;
            });
        }

        std::exception_ptr error = nullptr;

        for (auto &&pk : pkgs) {
            if (error) {
                // Wait for remaining futures
                pk.f_mechs.wait();
            } else {
                try {
                    pk.mechs = pk.f_mechs.get();
                } catch (...) {
                    error = std::current_exception();
                }
            }
        }

        if (error) {
            std::rethrow_exception(error);
        }

        packager.unpack(std::move(pkgs), geos, env);

        cat.write();
    }
}

int main(int argc, char *argv[]) {
    CommandLineArgs clargs = CommandLineArgs::parse(argc, argv);

    toml::v2::table config = toml::parse_file(clargs.config_file);

    std::unique_ptr<PotentialBase> ff = load_potential(config);

    auto [init, __] = load_supercell(config, ff->species_map());

    Classify classify = load_classifyer(config);

    Catalogue cat{options::Catalogue::load(config)};

    Packager packager(options::Packager::load(config));

    auto opt_find = options::FindMechanisms::load(config);

    double capt_tol = options::Mechanism::load(config).rel_cap_tol;

    std::vector<DiscreteKey> keys;
    std::vector<Geometry> geos;
    std::vector<Catalogue::pointer> env;

    std::unique_ptr minimise = load_minimiser(config);

    std::unique_ptr finder = load_sp_search(config);

    Visualise vis(options::Visualise::load(config));

    //
    auto m = tick("Minimise ");
    double Ei = ff->energy(init);
    bool initial_min = minimise->minimise(init, ff);
    double Ef = ff->energy(init);
    tock(m, "dE", Ei, Ef, Ef - Ei);

    ALWAYS_CHECK(initial_min, "Initial structure minimisation failed.");

    Supercell post_recon = init;

    Stream streamer(init);

    streamer.dump_raw(init, -1);

    update_catalogue(classify, cat, packager, init, finder, ff, opt_find, keys, geos, env);

    SuperCache superbasins{options::SuperCache::load(config), init, env};

    double time = 0;
    double time_lim = fetch<double>(config, "kinetics", "sim_time");

    VecN<double> ax = init.activ.zero_like();
    DimerRotor _rotor({.iter_max_rot = 200, .delta_r = 1e-3, .theta_tol = 1e-7});

    for (std::size_t i = 0; time < time_lim; i++) {
        /////////////////////////////////

        // Pre condition: caches active sb's occupied basin contains init

        auto [modified_cell, mech, dt, basin] = superbasins.select_mech(init);

        time += dt;

        if (modified_cell) {
            // Changed basin => Changed state => update geos
            // TODO: build only geo of mechs atom
            classify(init, keys, geos);
            streamer.dump_raw(init, -1);
        }

        /////////////////////////////////

        double E0 = ff->energy(init);

        auto const &m = superbasins.reconstruct(mech).onto(init, geos);

        post_recon.activ.view() = init.activ.view();

        double E1 = ff->energy(init);

        ALWAYS_CHECK(minimise->minimise(init, ff), "min_failed");

        /////////////////////////////////////////////////////////////

        try {
            update_catalogue(classify, cat, packager, init, finder, ff, opt_find, keys, geos, env);
        } catch (std::runtime_error const &error) {
            //
            std::cerr << error.what() << std::endl;

            streamer.dump_raw(init, -3);

            // Assuming init is in SP instead of min

            static thread_local pcg64 rng(pcg_extras::seed_seq_from<std::random_device>{});

            std::normal_distribution<> gauss_dist(0, 0.03);

            for (auto &&x : init.activ.view()) {
                x += gauss_dist(rng);
            }

            ALWAYS_CHECK(minimise->minimise(init, ff), "min_failed");

            streamer.dump_raw(init, -4);

            update_catalogue(classify, cat, packager, init, finder, ff, opt_find, keys, geos, env);
        }

        /////////////////////////////////

        double dR = norm(init.active_disp(post_recon.activ.view()));

        double Ef = ff->energy(init);

        std::cout << " It " << i                                                               //
                  << " ΔE " << m.activ_energy << "eV"                                          //
                  << " Ef " << Ef << "eV"                                                      //
                  << "  k " << m.pre_factor << "Hz"                                            //
                  << "  t " << dt << ":" << time << "s"                                        //
                  << " dE " << std::abs(E1 - Ef) << "eV"                                       //
                  << " dR " << dR << "A"                                                       //
                  << " dM " << m.rel_cap << ':' << m.abs_cap / m.rel_cap - m.abs_cap << '\n';  //

        streamer(init, i, time, E0, m.activ_energy, Ef, m.pre_factor);

        if (m.rel_cap <= capt_tol) {
            streamer.dump_raw(init, -2);
            std::cout << "Mech overflow, press ENTER to continue...\n";
            std::cin.ignore();
        }

        // Restore preconditions
        superbasins.connect_via(mech, init, env);

        std::cout << "\n";
    }

    cat.write();

    return 0;
}
