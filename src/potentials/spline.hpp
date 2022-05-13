#pragma once

#include <vector>

// Computes a set of cubic spline coefficients for a tabulated function, function/gradient then
// interpolated through appropriate methods
class Spline {
  public:
    // Interpolate tabulated function
    double operator()(double x) const;

    // Interpolate tabulated functions gradient
    double grad(double x) const;

    // Interpolate tabulated functions second derivative (not continuous/smooth)
    double grad2(double x) const;

    Spline() = default;

    // Construct from y, (n + 1) y_i values evenly spaced on interval 0,dx,...,ndx
    Spline(std::vector<double> y, double _dx);

  private:
    struct Spine {
        double a, b, c, d, cp, dp;
    };

    std::vector<Spine> _spines;

    double _dx;
    double _inv_dx;
};
