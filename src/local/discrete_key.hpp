#pragma once

#include <algorithm>
#include <array>
#include <vector>

#include "cereal/types/array.hpp"
#include "config.hpp"
#include "supercell.hpp"
#include "utility.hpp"

// Stores histograms of atom species & central atom. Introduces lexicographical ordering that can be
// used as a key to a std::map
struct DiscreteKey {
    Colour centre_col{};                   // Colour of central/first atom
    std::array<int, Colour::max()> sdf{};  // Species distribution function

    // Resets sdf to zero
    void clear() { std::fill(sdf.begin(), sdf.end(), 0); }

    // Strict weak ordering
    inline friend bool operator<(DiscreteKey const &a, DiscreteKey const &b) {
        if (a.centre_col != b.centre_col) {
            return a.centre_col < b.centre_col;
        }

        return std::lexicographical_compare(a.sdf.begin(), a.sdf.end(), b.sdf.begin(), b.sdf.end());
    }

    template <class Archive> void serialize(Archive &ar) { ar(centre_col, sdf); }
};
