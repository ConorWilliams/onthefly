#include "potentials/ADP/data.hpp"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <vector>

#include "potentials/spline.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

// Factory function, pareses tabulated adp data in LAMMPS adp/fs format
TabADP load_adp(toml::v2::table const &config) {
    auto in_file = fetch<std::string>(config, "potential", "in_file");

    std::ifstream file(in_file);

    ALWAYS_CHECK(file.good(), "Could not open adp file");

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

    TabADP data{cut};

    for (std::size_t i = 0; i < numS; i++) {
        std::string sp;
        line >> sp;

        auto [it, inserted] = data._map.emplace(std::move(sp), i);

        ALWAYS_CHECK(inserted, "Repeat species in EAM file");
    }

    std::vector<double> raw;

    for (std::size_t i = 0; i < numS; ++i) {
        // Read species info
        safe_getline(file) >> data.atomic(i) >> data.mass(i);

        // Read F
        for (std::size_t j = 0; j < numP; ++j) {
            double tmp;
            safe_getline(file) >> tmp;
            raw.push_back(tmp);
        }

        data.f(i) = Spline{raw, delP};
        raw.clear();

        // Read phi
        for (std::size_t j = 0; j < numR; ++j) {
            double tmp;
            safe_getline(file) >> tmp;
            raw.push_back(tmp);
        }

        data.phi(i) = Spline{raw, delR};
        raw.clear();
    }

    // Read v ***IMPORTANT**** tabulated as r*v
    for (std::size_t i = 0; i < numS; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            for (std::size_t k = 0; k < numR; ++k) {
                double tmp;
                safe_getline(file) >> tmp;
                double r = k * delR;
                raw.push_back(tmp / r);
            }

            raw[0] = 0;  // Fixup divide by zero for spline

            data.v(i, j) = Spline{raw, delR};
            raw.clear();
        }
    }

    // Read u(r): same format and symmetry as v
    for (std::size_t i = 0; i < numS; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            for (std::size_t k = 0; k < numR; ++k) {
                double tmp;
                safe_getline(file) >> tmp;
                raw.push_back(tmp);
            }

            raw[0] = 0;  // Fixup divide by zero for spline

            data.u(i, j) = Spline{raw, delR};
            raw.clear();
        }
    }

    // Read w(r): same format and symmetry as v
    for (std::size_t i = 0; i < numS; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            for (std::size_t k = 0; k < numR; ++k) {
                double tmp;
                safe_getline(file) >> tmp;
                raw.push_back(tmp);
            }

            raw[0] = 0;  // Fixup divide by zero for spline

            data.w(i, j) = Spline{raw, delR};
            raw.clear();
        }
    }

    // Test
    // std::ofstream testfile;
    // testfile.open("data.txt");
    // testfile << "numS: " << numS << "\n";
    // testfile << "numP: " << numP << ", delP: " << delP << ", numR: " << numR << ", delR: " <<
    // delR
    //          << ", Cut: " << cut << "\n";
    // testfile << "F(rho) Cr \t\t Rho(r) Cr \t\t F(Rho) Ni \t\t Rho(r) Ni \t\t V(r) Cr Cr \t\t V(r)
    // "
    //             "Ni Cr \t\t V(r) Ni Ni \t\t u(r) Cr Cr \t\t u(r) Ni Cr \t\t u(r) Ni Ni \t\t w(r)
    //             " "Cr Cr \t\t w(r) Ni Cr \t\t w(r) Ni Ni \n";
    // for (std::size_t i = 0; i < numR; i++) {
    //     testfile << data.f(0)(i * delP) << "\t\t" << data.phi(0)(i * delR) << "\t\t"
    //              << data.f(1)(i * delP) << "\t\t" << data.phi(1)(i * delR) << "\t\t"
    //              << data.v(0, 0)(i * delR) << "\t\t" << data.v(1, 0)(i * delR) << "\t\t"
    //              << data.v(1, 1)(i * delR) << "\t\t" << data.u(0, 0)(i * delR) << "\t\t"
    //              << data.u(1, 0)(i * delR) << "\t\t" << data.u(1, 1)(i * delR) << "\t\t"
    //              << data.w(0, 0)(i * delR) << "\t\t" << data.w(1, 0)(i * delR) << "\t\t"
    //              << data.w(1, 1)(i * delR) << "\n";
    // }
    // testfile.close();
    // End Test

    return data;
}
