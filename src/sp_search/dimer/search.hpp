#pragma once

#include <memory>
#include <random>

#include "minimise/minimiser_base.hpp"
#include "pcg_random.hpp"
#include "potentials/potential_base.hpp"
#include "sp_search/dimer/dimer.hpp"
#include "sp_search/sp_search_base.hpp"
#include "sp_search/vineyard.hpp"
#include "supercell.hpp"

template <typename T> class DimerSPS final : public SearchBase {
  public:
    DimerSPS(std::unique_ptr<MinimiserBase> minimiser, T const& dimer);

    virtual std::unique_ptr<SearchBase> clone() const override;

    bool find_sp(Supercell const& init,
                 Supercell& dimer,
                 Supercell& final,
                 std::unique_ptr<PotentialBase>& ff) override;

  private:
    T _dimer;
    std::unique_ptr<MinimiserBase> _minimiser;

    double _nudge;
    double _tol;

    VecN<double> _ax;

    Supercell _old;
    // DimerRotor _rotor{{6, 1000, 0.01, 0.0001}};
};

template <typename T>
DimerSPS<T>::DimerSPS(std::unique_ptr<MinimiserBase> minimiser, T const& dimer)
    : _dimer(dimer),
      _minimiser(std::move(minimiser)),
      _nudge(dimer.get_opt().nudge),
      _tol(dimer.get_opt().basin_tol) {}

template <typename T> std::unique_ptr<SearchBase> DimerSPS<T>::clone() const {
    return std::make_unique<DimerSPS<T>>(_minimiser->clone(), _dimer);
}

template <typename T> bool DimerSPS<T>::find_sp(Supercell const& init,
                                                Supercell& dimer,
                                                Supercell& final,
                                                std::unique_ptr<PotentialBase>& ff) {
    // Build dimer axis
    _ax = dimer.activ.view() - init.activ.view();

    // Randomises the elements in v (multiplicatively) then normalises.
    static thread_local pcg64 rng(pcg_extras::seed_seq_from<std::random_device>{});

    std::normal_distribution<> gauss_dist(0, 1);

    for (int i = 0; i < _ax.size(); ++i) {
        _ax[i] *= gauss_dist(rng);
    }

    _ax *= 1 / norm(_ax);

    //////////////////////////////////////////////

    // Do SPS
    if (!_dimer.find_sp(dimer, _ax, ff)) {
        return false;
    }

    final = dimer;
    _old = dimer;

    final.activ.view() += _ax * _nudge;
    _old.activ.view() -= _ax * _nudge;

    if (!_minimiser->minimise(_old, ff, _ax) || !_minimiser->minimise(final, ff, _ax)) {
        std::cout << "Minimisation failed\n";
        return false;
    }

    double final_init = norm(final.active_disp(init.activ.view()));
    double old_init = norm(_old.active_disp(init.activ.view()));

    if (old_init > final_init) {
        using std::swap;
        swap(old_init, final_init);
        swap(final, _old);
    }

    if (old_init > _tol) {
        // std::cout << "F1 old-init (disconnected): " << old_init << '\n';
        return false;
    }

    if (norm(_old.active_disp(final.activ.view())) < _tol) {
        // std::cout << "F2 old-fin: " << norm(_old.active_disp(final.activ.view())) << '\n';
        return false;
    }

    if (norm(_old.active_disp(dimer.activ.view())) < _tol * 0.25) {
        std::cout << "F3 old-dim: " << norm(_old.active_disp(dimer.activ.view())) << '\n';
        return false;
    }

    if (norm(final.active_disp(dimer.activ.view())) < _tol * 0.25) {
        std::cout << "F4 fin-dim: " << norm(final.active_disp(dimer.activ.view())) << '\n';
        return false;
    }

    // if (_rotor.align(final, _ax, ff) < 0) {
    //     // std::cout << "Converged to sp\n";
    //     return false;
    // }

    // dump_supercell(_old, "olkmc.xyz", true);
    // dump_supercell(dimer, "olkmc.xyz", true);
    // dump_supercell(final, "olkmc.xyz", true);
    // dump_supercell(final, "fin.xyz", true);

    // exit(1);

    return true;
}
