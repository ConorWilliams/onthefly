#pragma once

#include "config.hpp"
#include "potentials/potential_base.hpp"
#include "supercell.hpp"

class Vineyard {
  public:
    Vineyard(double zero_tol) : _tol(zero_tol) {}

    // Compute the Eigen-values of the hessian at a minima
    void load_basin(Supercell const&, std::unique_ptr<PotentialBase>&);

    // Compute the Eigen-values of the hessian at a sp, returns false if not 1st order
    bool load_sp(Supercell const&, std::unique_ptr<PotentialBase>&);

    // Compute the harmonic prefactor of the loaded basin/sp supercells
    double pre_factor();

  private:
    double _tol;

    VecN<double> _ev_basin;
    VecN<double> _ev_sp;
};
