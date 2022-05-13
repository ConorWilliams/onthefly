#include "kinetics/supercache.hpp"

#include <bits/stdint-intn.h>

#include <iostream>
#include <random>

#include "kinetics/basin.hpp"

namespace options {

SuperCache SuperCache::load(toml::v2::table const &config) {
    SuperCache opt;

    opt.state_tol = fetch<double>(config, "kinetics", "state_tol");
    opt.barrier_tol = fetch<double>(config, "kinetics", "barrier_tol");
    opt.cache_size = config["kinetics"]["cache_size"].value_or(opt.cache_size);

    opt.dynamic_tol = config["kinetics"]["dynamic_tol"].value_or(opt.dynamic_tol);

    if (opt.dynamic_tol) {
        opt.max_superbasin_size = fetch<int64_t>(config, "kinetics", "max_superbasin_size");
        opt.tol_grow = fetch<double>(config, "kinetics", "tol_grow");
        opt.tol_shrink = fetch<double>(config, "kinetics", "tol_shrink");
    }

    opt.basin = Basin::load(config);

    return opt;
}

}  // namespace options

void SuperCache::connect_via(std::size_t mech,
                             Supercell const &cell,
                             std::vector<Catalogue::pointer> const &env) {
    // std::cout << "low barrier" << (*_sb)[mech].barrier << ' ' << _opt.barrier_tol << "\n";
    if (std::optional<std::size_t> basin = _sb.find_occupy(cell, _opt.state_tol)) {
        // Try and connect states if did internal jump
        _sb.connect_from(*basin, mech);
        std::cout << "Existing basin in SB, size: " << _sb.size() << std::endl;

    } else if ((*_sb)[mech].barrier < _opt.barrier_tol) {
        // Followed low barrier to get here
        if (_opt.dynamic_tol && _sb.size() >= _opt.max_superbasin_size) {
            // Overflowing SB, must lower barrier_tol
            _opt.barrier_tol = std::max(0.0, _opt.barrier_tol * _opt.tol_shrink);
            _sb = Superbasin({_opt.basin, cell, env});
            _cache.clear();

            std::cout << "Dynamically adjusting, barrier_tol=" << _opt.barrier_tol << '\n';
        } else {
            // Can expand and occupy;
            _sb.connect_from(_sb.expand_occupy({_opt.basin, cell, env}), mech);
            std::cout << "New basin in SB, size: " << _sb.size() << std::endl;
        }
    } else {
        // Therefore followed high-barrier out of basin

        // Try and retrieve cached SB
        std::optional cached = [&]() -> std::optional<Superbasin> {
            for (auto it = _cache.begin(); it != _cache.end(); ++it) {
                if (std::optional basin = it->find_occupy(cell, _opt.state_tol)) {
                    Superbasin tmp = std::move(*it);
                    _cache.erase(it);
                    return tmp;
                }
            }
            return std::nullopt;
        }();

        if (cached) {
            std::cout << "!======LOAD_CACHED=====! " << cached->size() << ':' << size() << "\n";
            cache(std::exchange(_sb, std::move(*cached)));
            _in_cache_count += 1;
        } else {
            std::cout << "!=======NEW_SUPER======! " << size() << "\n";
            cache(std::exchange(_sb, Superbasin({_opt.basin, cell, env})));
            _in_cache_count = 0;
        }

        if (_opt.dynamic_tol && _in_cache_count > _opt.cache_size) {
            _opt.barrier_tol *= _opt.tol_grow;
            _sb = Superbasin({_opt.basin, cell, env});
            _cache.clear();

            std::cout << "Dynamically adjusting, barrier_tol=" << _opt.barrier_tol << '\n';
        }
    }
}

namespace {

thread_local pcg64 psudo_rng{pcg_extras::seed_seq_from<std::random_device>{}};

}  // namespace

Basin::choice SuperCache::select_mech(Supercell &cell) {
    if ((*_sb).connected) {
        Basin::choice choice = _sb.kmc_choice(psudo_rng);

        if (choice.basin_changed) {
            cell.activ.view() = (*_sb).state();
        }

        return choice;
    } else {
        return (*_sb).kmc_choice(psudo_rng, _sb.occupied());
    }
}

// If space, enough cache SB and delete old SBs if necessary
void SuperCache::cache(Superbasin &&basin) {
    //
    _cache.push_front(std::move(basin));

    if (_cache.size() > _opt.cache_size) {
        _cache.pop_back();
    }
}
