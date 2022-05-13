#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>
#include <vector>

#include "cereal/types/array.hpp"
#include "cereal/types/vector.hpp"
#include "config.hpp"
#include "supercell.hpp"
#include "utility.hpp"

class fuzzy_key {
  public:
    template <class Atom_t> void build(std::vector<Atom_t> const& atoms) {
        for (std::size_t i = 1; i < atoms.size(); i++) {
            _cen[atoms[i].col].push_back(norm(atoms[i].vec - atoms[0].vec));

            for (std::size_t j = 1; j < i; j++) {
                _fuzzy(atoms[i].col, atoms[j].col).push_back(norm(atoms[i].vec - atoms[j].vec));
            }
        }

        // Put _fuzzy_key into sorted order
        for (auto&& vec : _cen) {
            std::sort(vec.begin(), vec.end());
        }
        _fuzzy.foreach ([](auto& vec) { std::sort(vec.begin(), vec.end()); });
    }

    void clear() {
        for (auto&& vec : _cen) {
            vec.clear();
        }
        _fuzzy.foreach ([](auto& vec) { vec.clear(); });
    }

    // If keys a and b equivalent within delta returns true
    friend bool equivalent(double delta, fuzzy_key const& a, fuzzy_key const& b) {
        //

        for (std::size_t i = 0; i < Colour::max(); i++) {
            if (a._cen[i].size() != b._cen[i].size()) {
                return false;
            }

            for (std::size_t k = 0; k < a._cen[i].size(); ++k) {
                if (std::abs(a._cen[i][k] - b._cen[i][k]) > SQRT_2 * delta) {
                    return false;
                }
            }

            for (std::size_t j = 0; j <= i; j++) {
                if (a._fuzzy(i, j).size() != b._fuzzy(i, j).size()) {
                    return false;
                }

                for (std::size_t k = 0; k < a._fuzzy(i, j).size(); ++k) {
                    if (std::abs(a._fuzzy(i, j)[k] - b._fuzzy(i, j)[k]) > SQRT_2 * delta) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    // Returns L-inf norm, between fuzzy keys
    friend double fuzzy_norm(fuzzy_key const& a, fuzzy_key const& b) {
        //
        double max = 0;

        for (std::size_t i = 0; i < Colour::max(); i++) {
            CHECK(a._cen[i].size() == b._cen[i].size(), "fkeys wrong size A");

            for (std::size_t k = 0; k < a._cen[i].size(); ++k) {
                max = std::max(max, std::abs(a._cen[i][k] - b._cen[i][k]));
            }

            for (std::size_t j = 0; j <= i; j++) {
                CHECK(a._fuzzy(i, j).size() == b._fuzzy(i, j).size(), "fkeys wrong size B");

                for (std::size_t k = 0; k < a._fuzzy(i, j).size(); ++k) {
                    max = std::max(max, std::abs(a._fuzzy(i, j)[k] - b._fuzzy(i, j)[k]));
                }
            }
        }

        return max;
    }

    template <class Archive> void serialize(Archive& ar) { ar(_cen, _fuzzy); }

  private:
    std::array<std::vector<double>, Colour::max()> _cen{};
    SymMat<std::vector<double>, Colour::max()> _fuzzy{};
};

namespace experimental {

// Experimental (overfuzzed) key
class fuzzy_key {
  public:
    template <class Atom_t> void build(std::vector<Atom_t> const& atoms) {
        // Initialise
        _fuzzy.resize(atoms.size(), std::vector<double>(atoms.size()));

        for (std::size_t i = 0; i < atoms.size(); i++) {
            for (std::size_t j = 0; j < atoms.size(); j++) {
                _fuzzy[i][j] = norm(atoms[i].vec - atoms[j].vec);
            }

            std::sort(_fuzzy[i].begin(), _fuzzy[i].end());
        }

        // Put _fuzzy_key into sorted order
        std::sort(_fuzzy.begin() + 1, _fuzzy.end(), [](auto const& a, auto const& b) {
            // return a[0] < b[0];
            return std::reduce(a.begin(), a.end()) < std::reduce(b.begin(), b.end());
        });
    }

    void clear() {
        for (auto&& vec : _fuzzy) {
            vec.clear();
        }
    }

    inline friend bool equivalent(double delta, fuzzy_key const& a, fuzzy_key const& b) {
        if (a._fuzzy.size() != b._fuzzy.size()) {
            return false;
        }

        for (std::size_t i = 0; i < a._fuzzy.size(); i++) {
            for (std::size_t j = 0; j < a._fuzzy[i].size(); j++) {
                if (std::abs(a._fuzzy[i][j] - b._fuzzy[i][j]) > SQRT_2 * delta) {
                    return false;
                }
            }
        }

        return true;
    }

    template <class Archive> void serialize(Archive& ar) { ar(_fuzzy); }

  private:
    std::vector<std::vector<double>> _fuzzy;
};

}  // namespace experimental
