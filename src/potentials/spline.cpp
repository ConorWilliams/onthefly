#include "potentials/spline.hpp"

#include <cstddef>

#include "utility.hpp"

// Computes a set of cubic spline coefficients for a tabulated function, function/gradient then
// interpolated through appropriate methods
double Spline::operator()(double x) const {
    CHECK(x >= 0, "x out of bounds");

    std::size_t const i = x * _inv_dx;

    CHECK(i < _spines.size(), "i out of bounds");

    x -= i * _dx;

    return _spines[i].a + x * (_spines[i].b + x * (_spines[i].c + x * _spines[i].d));
}

// Interpolate tabulated functions gradient
double Spline::grad(double x) const {
    CHECK(x >= 0, "x out of bounds ");

    std::size_t const i = x * _inv_dx;

    CHECK(i < _spines.size(), "i out of bounds in grad");

    x -= i * _dx;

    return _spines[i].b + x * (_spines[i].cp + x * _spines[i].dp);
}

// Interpolate tabulated functions gradient
double Spline::grad2(double x) const {
    CHECK(x >= 0, "x out of bounds ");

    std::size_t const i = x * _inv_dx;

    CHECK(i < _spines.size(), "i out of bounds in grad");

    x -= i * _dx;

    return _spines[i].cp + 2.0 * x * _spines[i].dp;
}

// Based on Wikipedia algorithm: https://en.wikipedia.org/wiki/Spline_(mathematics)
// y is (n + 1) y_i values evenly spaced on interval 0,dx,...,ndx
Spline::Spline(std::vector<double> y, double _dx) : _spines{}, _dx(_dx), _inv_dx(1 / _dx) {
    //

    y.push_back(double(y.back()));  // pad spline

    std::size_t n = y.size() - 1;

    // 1
    std::vector<double> a(n + 1);
    // 2
    std::vector<double> b(n);
    std::vector<double> d(n);
    // 4
    std::vector<double> alpha(n);
    // 5
    std::vector<double> c(n + 1);
    std::vector<double> l(n + 1);
    std::vector<double> mu(n + 1);
    std::vector<double> z(n + 1);

    // 1
    for (std::size_t i = 0; i <= n; ++i) {
        a[i] = y[i];
    }

    // 3
    // h_i = _dx

    // 4
    for (std::size_t i = 1; i <= n - 1; ++i) {
        const double inv_h = 3 / _dx;

        alpha[i] = inv_h * (a[i + 1] - 2 * a[i] + a[i - 1]);
    }

    // 6
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    // 7
    for (std::size_t i = 1; i <= n - 1; ++i) {
        l[i] = _dx * (4 - mu[i - 1]);
        mu[i] = _dx / l[i];
        z[i] = (alpha[i] - _dx * z[i - 1]) / l[i];
    }

    // 8
    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    // 9
    for (std::size_t i = 1; i <= n; ++i) {
        std::size_t j = n - i;

        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / _dx - _dx * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * _dx);
    }

    // 11
    for (std::size_t i = 0; i <= n - 1; ++i) {
        _spines.push_back({a[i], b[i], c[i], d[i], 2 * c[i], 3 * d[i]});
    }
}
