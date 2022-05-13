

#include "potentials/EAM/data.hpp"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

#include "potentials/spline.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

// Read N elements (doubles) in lines of length K
std::vector<double> read_chunked(std::ifstream& file, std::size_t n, std::size_t k) {
    //
    std::vector<double> raw;
    std::istringstream line;

    for (std::size_t i = 0; i < n; ++i) {
        if (i % k == 0) {
            line = safe_getline(file);
        }

        double tmp;
        line >> tmp;
        raw.push_back(tmp);
    }

    return raw;
}

// Factory function, pareses tabulated eam data in LAMMPS eam/fs format
TabEAM load_eam(toml::v2::table const& config) {
    //
    std::ifstream file(fetch<std::string>(config, "potential", "in_file"));

    ALWAYS_CHECK(file.good(), "Could not open eam file");

    // Temporaries
    std::size_t numS, numP, numR;
    double delP, delR, cut;

    // Skip to 4th line
    for (int i = 0; i < 3; ++i) {
        safe_getline(file);
    }

    // Get number of species
    auto line = safe_getline(file);

    line >> numS;

    ALWAYS_CHECK(numS == NUM_ATOM_SPECIES, "Compile time number of species != dynamic count");

    safe_getline(file) >> numP >> delP >> numR >> delR >> cut;

    TabEAM data{cut};

    for (std::size_t i = 0; i < numS; i++) {
        std::string sp;
        line >> sp;

        auto [it, inserted] = data._map.emplace(std::move(sp), i);

        ALWAYS_CHECK(inserted, "Repeat species in EAM file");
    }

    for (std::size_t i = 0; i < numS; ++i) {
        // Read species info
        safe_getline(file) >> data.atomic(i) >> data.mass(i);

        // Read F
        data.f(i) = Spline{read_chunked(file, numP, 5), delP};

        // Read phi
        for (std::size_t j = 0; j < numS; ++j) {
            data.phi(i, j) = Spline{read_chunked(file, numP, 5), delR};
        }
    }

    // Read v ***IMPORTANT**** tabulated as r*v
    for (std::size_t i = 0; i < numS; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            auto raw = read_chunked(file, numP, 5);

            for (size_t i = 0; i < raw.size(); i++) {
                raw[i] /= delR * i;
            }

            ALWAYS_CHECK(raw.size() > 0, "no elements!");

            raw[0] = 0;  // Fixup divide by zero for spline

            data.v(i, j) = Spline{std::move(raw), delR};
        }
    }

    return data;
}