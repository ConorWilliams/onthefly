#include "kinetics/superbasin.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <utility>

#include "config.hpp"
#include "pcg_random.hpp"

Basin::choice Superbasin::kmc_choice(pcg64 &psudo_rng) {
    //
    std::uniform_real_distribution<double> uniform;

    VecN<double> tau = compute_tau();

    int count = 0;

    // Sum over all basin->escape rate times basin modifiers, omit normalising factor of 1/tau
    double const r_sum = [&]() {
        double sum = 0;

        for (std::size_t i = 0; i < size(); ++i) {
            double basin_exit_sum = 0;
            for (std::size_t j = 0; j < _super[i].size(); ++j) {
                // basin -> non-basin
                if (_super[i][j].exit_mech) {
                    ++count;
                    basin_exit_sum += _super[i][j].rate;
                }
            }
            sum += tau[i] * basin_exit_sum;
        }
        return sum;
    }();

    CHECK(count > 0, "o no countin wong");
    CHECK(r_sum > 0, "r_sum negative");

    auto [basin, mech] = [&]() -> std::pair<std::size_t, std::size_t> {
        double const lim = uniform(psudo_rng) * r_sum;
        double sum = 0;

        for (std::size_t i = 0; i < size(); ++i) {
            for (std::size_t j = 0; j < _super[i].size(); ++j) {
                if (_super[i][j].exit_mech) {
                    sum += tau[i] * _super[i][j].rate;
                    if (sum > lim) {
                        return {i, j};
                    }
                }
            }
        }
        throw std::runtime_error("Hit end of super choice");
    }();

    // Must occupy new basin (this is why method is not const)
    std::size_t old_basin = std::exchange(_occupied, basin);

    double const inv_tau = 1 / tau.sum();
    double const eff_rate = tau[basin] * inv_tau * _super[basin][mech].rate;
    double const prob = 100 * eff_rate / (inv_tau * r_sum);

    std::cout << "Effective rate: " << eff_rate << " : " << prob << "% of " << count
              << "  choices\n";

    // Must normalize by inv_tau
    return {
        old_basin != basin,
        mech,
        -std::log(uniform(psudo_rng)) / (r_sum * inv_tau),
        _occupied,
    };
}

// Compute vector of mean residence-times in each basin
VecN<double> Superbasin::compute_tau() const {
    // Non-allocating Eigen3-objects: theta_{i} = Kroneker_{is}
    auto lam = [s = _occupied](std::size_t i) { return i == s; };
    auto theta = VecN<double>::NullaryExpr(size(), std::move(lam));
    auto identity = MatN<double>::Identity(size(), size());

    // Calc theta^{sum}
    VecN<double> tau = (identity - _prob).colPivHouseholderQr().solve(theta.matrix());

    // Convert to tau_i
    for (std::size_t i = 0; i < size(); ++i) {
        tau[i] /= _super[i].rate_sum();
    }

    return tau;
}

void Superbasin::connect_from(std::size_t i, std::size_t m) {
    _prob(_occupied, i) = _super[i][m].rate / _super[i].rate_sum();
    _super[i][m].exit_mech = false;
    _super[i].connected = true;
    return;
}

std::size_t Superbasin::expand_occupy(Basin &&basin) {
    _super.push_back(std::move(basin));
    _prob.conservativeResizeLike(MatN<double>::Zero(size(), size()));
    return std::exchange(_occupied, size() - 1);
}

std::optional<std::size_t> Superbasin::find_occupy(Supercell const &cell, double tol) {
    for (std::size_t i = 0; i < size(); ++i) {
        if (norm_sq(cell.active_disp(_super[i].state())) < tol * tol) {
            return std::exchange(_occupied, i);
        }
    }
    return std::nullopt;
}
