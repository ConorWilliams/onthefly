#include "sp_search/vineyard.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>

#include "Eigen/Eigenvalues"
#include "minimise/LBFGS/lbfgs.hpp"
#include "supercell.hpp"

namespace {

void compute_eigen_values(Supercell const& cell,
                          std::unique_ptr<PotentialBase>& ff,
                          VecN<double>& out) {
    //

    static thread_local MatN<double> hess;
    static thread_local Eigen::SelfAdjointEigenSolver<MatN<double>> es;

    ff->hessian(cell, hess);

    es.compute(hess, Eigen::EigenvaluesOnly);

    out = es.eigenvalues();
}

}  // namespace

// Compute the Eigen-values of the hessian at a minima
void Vineyard::load_basin(Supercell const& cell, std::unique_ptr<PotentialBase>& ff) {
    //

    // dump_supercell(cell, "vine_fail.xyz", true);

    compute_eigen_values(cell, ff, _ev_basin);

    if (auto order = (_ev_basin < -_tol).count(); order != 0) {
        //
        std::cout << "** WARN ** " << _ev_basin.transpose().head(10) << std::endl;

        throw std::runtime_error("Vineyard load_basin contract violated");

        // static thread_local VecN<double> out;

        // ff->gradient(cell, out);

        // std::cout << "FORCE = " << norm(out) << std::endl;

        // ALWAYS_CHECK(false, "Not a basin: " + std::to_string(order));
    }
}

// Compute the Eigen-values of the hessian at a sp
bool Vineyard::load_sp(Supercell const& cell, std::unique_ptr<PotentialBase>& ff) {
    //
    compute_eigen_values(cell, ff, _ev_sp);

    std::size_t order = (_ev_sp < -_tol).count();

    switch (order) {
        case 0:
            std::cout << _ev_sp.transpose().head(10) << std::endl;
            std::cout << "SP is minima\n";
            return false;
        case 1:
            return true;
        default:
            std::cout << _ev_sp.transpose().head(10) << std::endl;
            std::cout << "SP is " << order << " order\n";
            return false;
    }
}

// Compute the harmonic prefactor of the loaded basin/sp supercells
double Vineyard::pre_factor() {
    //
    CHECK(_ev_sp.size() > 0 && _ev_sp.size() == _ev_basin.size(), "Incorrectly primed Vine");

    double prod = 1;

    for (int i = 0; i < _ev_sp.size(); i++) {
        if (_ev_basin[i] > _tol) {
            prod *= _ev_basin[i];
        }
        if (_ev_sp[i] > _tol) {
            prod /= _ev_sp[i];
        }
    }
    // Correcting for mass in AMU: 1

    prod = std::sqrt(prod / (2 * M_PI * 1.6605390666050e-27));

    // std::cout << "The harmonic pre-factor is: " << prod << "Hz\n";

    return prod;
}
