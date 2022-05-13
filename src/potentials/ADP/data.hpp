#pragma once

#include <array>
#include <cstddef>
#include <unordered_map>

#include "config.hpp"
#include "potentials/spline.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

// Holds tabulated TabADP data for N atom species
class TabADP {
  public:
    double rcut() const { return _rcut; }
    double rcut() { return _rcut; }

    Spline const &f(std::size_t i) const {
        CHECK(i < N, "Bad Access");
        return _f[i];
    }
    Spline const &phi(std::size_t i) const {
        CHECK(i < N, "Bad Access");
        return _phi[i];
    }
    Spline const &v(std::size_t i, std::size_t j) const {
        CHECK(i < N && j < N, "Bad Access");
        return _v(i, j);
    }
    Spline const &u(std::size_t i, std::size_t j) const {
        CHECK(i < N && j < N, "Bad Access");
        return _u(i, j);
    }
    Spline const &w(std::size_t i, std::size_t j) const {
        CHECK(i < N && j < N, "Bad Access");
        return _w(i, j);
    }
    double const &mass(std::size_t i) const {
        CHECK(i < N, "Bad Access");
        return _mass[i];
    }
    std::size_t const &atomic(std::size_t i) const {
        CHECK(i < N, "Bad Access");
        return _atomic[i];
    }

    std::unordered_map<std::string, std::uint16_t> const &species_map() const { return _map; }

  private:
    static constexpr std::size_t N = NUM_ATOM_SPECIES;

    double _rcut;

    std::array<double, N> _mass;
    std::array<std::size_t, N> _atomic;

    std::array<Spline, N> _f;
    std::array<Spline, N> _phi;
    SymMat<Spline, N> _v;
    SymMat<Spline, N> _u;
    SymMat<Spline, N> _w;

    std::unordered_map<std::string, std::uint16_t> _map;

    friend TabADP load_adp(toml::v2::table const &config);

    TabADP(double rcut) : _rcut(rcut) {}

    Spline &f(std::size_t i) {
        CHECK(i < N, "Bad Access");
        return _f[i];
    }
    Spline &phi(std::size_t i) {
        CHECK(i < N, "Bad Access");
        return _phi[i];
    }
    Spline &v(std::size_t i, std::size_t j) {
        CHECK(i < N && j < N, "Bad Access");
        return _v(i, j);
    }
    Spline &u(std::size_t i, std::size_t j) {
        CHECK(i < N && j < N, "Bad Access");
        return _u(i, j);
    }
    Spline &w(std::size_t i, std::size_t j) {
        CHECK(i < N && j < N, "Bad Access");
        return _w(i, j);
    }
    double &mass(std::size_t i) {
        CHECK(i < N, "Bad Access");
        return _mass[i];
    }
    std::size_t &atomic(std::size_t i) {
        CHECK(i < N, "Bad Access");
        return _atomic[i];
    }
};

// Factory function, pareses tabulated adp data in LAMMPS adp/fs format
TabADP load_adp(toml::v2::table const &config);
