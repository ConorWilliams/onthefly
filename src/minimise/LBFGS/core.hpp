#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "config.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

// LBFGS abstraction, holds variable history and computes the next newton step at each call
class CoreLBFGS {
  public:
    CoreLBFGS(std::size_t n) : _n{n}, _k{0}, _s(n), _y(n), _rho(n), _a(n){};

    // Reset core for re-use
    void clear() { _k = 0; }

    // Computes the q = H_k g_k using the L-BFGS two-loop recursion
    template <class E> void operator()(E const &x, VecN<double> const &gx, VecN<double> &q);

  private:
    std::size_t _n;
    std::size_t _k;

    std::vector<VecN<double>> _s;
    std::vector<VecN<double>> _y;

    VecN<double> _rho;
    VecN<double> _a;

    VecN<double> _prev_x;
    VecN<double> _prev_gx;
};

// Computes the q = H_k g_k using the L-BFGS two-loop recursion
template <class E> void CoreLBFGS::operator()(E const &x, VecN<double> const &gx, VecN<double> &q) {
    std::size_t prev = (_k - 1) % _n;

    // Compute the k-1 th  y, s and rho
    if (_k > 0) {
        using std::swap;

        swap(_s[prev], _prev_x);
        _s[prev] -= x;

        swap(_y[prev], _prev_gx);
        _y[prev] -= gx;

        // If Wolfie conditions then dot(y, s) > 0, if bad line search and wolfie not fullfilled we
        // take absolute value to prevent ascent direction

        _rho[prev] = 1.0 / std::abs(dot(_s[prev], _y[prev]));
    }

    _prev_x = x;
    _prev_gx = gx;

    int incur = _k <= _n ? 0 : _k - _n;
    int bound = _k <= _n ? _k : _n;

    q = gx;

    // Loop 1
    for (int i = bound - 1; i >= 0; --i) {
        int j = (i + incur) % _n;

        _a[j] = _rho[j] * dot(_s[j], q);
        q -= _a[j] * _y[j];
    }

    // Scaling Hessian_0
    if (_k > 0) {
        q *= 1 / (_rho[prev] * dot(_y[prev], _y[prev]));
    }

    // Loop 2
    for (int i = 0; i <= bound - 1; ++i) {
        int j = (i + incur) % _n;

        double b = _rho[j] * dot(_y[j], q);

        q += (_a[j] - b) * _s[j];
    }

    ++_k;
}
