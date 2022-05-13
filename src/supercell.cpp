#include "supercell.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>

#include "config.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

Simbox::Simbox(Vec3<double> const &extents, std::size_t periodic)
    : extents(extents), inv_extents(1.0 / extents), periodic(periodic) {
    ALWAYS_CHECK((extents > 0).all(), "Simbox extents are negative");
}

// Returns the smallest vector connecting a to a periodic image of b
Vec3<double> Simbox::min_image(Vec3<double> const &a, Vec3<double> const &b) const {
    using namespace periodic;

    Vec3<double> dr = b - a;

    switch (periodic) {
        case x | y | z:
            return dr - extents * (dr * inv_extents + 0.5).floor();
        case x | y:
            return {
                dr[0] - extents[0] * std::floor(dr[0] * inv_extents[0] + 0.5),
                dr[1] - extents[1] * std::floor(dr[1] * inv_extents[1] + 0.5),
                dr[2],
            };
        case y | z:
            return {
                dr[0],
                dr[1] - extents[1] * std::floor(dr[1] * inv_extents[1] + 0.5),
                dr[2] - extents[2] * std::floor(dr[2] * inv_extents[2] + 0.5),
            };
        case x | z:
            return {
                dr[0] - extents[0] * std::floor(dr[0] * inv_extents[0] + 0.5),
                dr[1],
                dr[2] - extents[2] * std::floor(dr[2] * inv_extents[2] + 0.5),
            };
        case x:
            return {
                dr[0] - extents[0] * std::floor(dr[0] * inv_extents[0] + 0.5),
                dr[1],
                dr[2],
            };
        case y:
            return {
                dr[0],
                dr[1] - extents[1] * std::floor(dr[1] * inv_extents[1] + 0.5),
                dr[2],
            };
        case z:
            return {
                dr[0],
                dr[1],
                dr[2] - extents[2] * std::floor(dr[2] * inv_extents[2] + 0.5),
            };
        case non:
            return dr;
        default:
            ALWAYS_CHECK(false, "Invalid periodicity");
    }
}

Simbox Simbox::load(toml::v2::table const &config) {
    auto lx = fetch<double>(config, "supercell", "simbox", "lx");
    auto ly = fetch<double>(config, "supercell", "simbox", "ly");
    auto lz = fetch<double>(config, "supercell", "simbox", "lz");

    std::size_t p = periodic::non;

    if (fetch<bool>(config, "supercell", "simbox", "px")) {
        p |= periodic::x;
    }
    if (fetch<bool>(config, "supercell", "simbox", "py")) {
        p |= periodic::y;
    }
    if (fetch<bool>(config, "supercell", "simbox", "pz")) {
        p |= periodic::z;
    }

    return {{lx, ly, lz}, p};
}

VecN<double> Supercell::active_disp(VecN<double> const &other) const {
    CHECK(activ.view().size() == other.size(), "Number of atoms are different");

    VecN<double> result = other - activ.view();

    // Can skip translational correction if bound atoms are present as they prevent drift
    if (bound.size() > 0) {
        return result;
    }

    // Find COM of atoms close together in both sets. If they are different ends of a mechanism with
    // a translational component this should give preference to atoms NOT involved in mechanism

    std::size_t count = 0;

    Vec3<double> com_a = Vec3<double>::Zero();
    Vec3<double> com_b = Vec3<double>::Zero();

    for (int i = 0; i < activ.view().size() / 3; ++i) {
        com_a[0] += activ.view()[3 * i + 0];
        com_a[1] += activ.view()[3 * i + 1];
        com_a[2] += activ.view()[3 * i + 2];

        com_b[0] += other[3 * i + 0];
        com_b[1] += other[3 * i + 1];
        com_b[2] += other[3 * i + 2];

        count += 1;
    }

    // std::cout << min << std::endl;

    ALWAYS_CHECK(count > 0, "divide by zero");

    com_a /= count;
    com_b /= count;

    // Remove translational component

    Vec3<double> delta = com_b - com_a;

    for (int i = 0; i < activ.view().size() / 3; ++i) {
        result[3 * i + 0] -= delta[0];
        result[3 * i + 1] -= delta[1];
        result[3 * i + 2] -= delta[2];
    }

    return result;
}

// In LAMMPS compatable .xyz file
void dump_supercell(Supercell const &cell, std::string const &out_file, bool append) {
    //
    std::ofstream outfile;
    if (append) {
        outfile.open(out_file, std::ios_base::app);
    } else {
        outfile.open(out_file);
    }

    outfile << cell.size() << std::endl;
    outfile << "Lattice=\"" << cell.extents[0] << " 0.0 0.0 0.0 " << cell.extents[1]
            << " 0.0 0.0 0.0 " << cell.extents[2] << "\"";

    outfile << std::setprecision(15);

    //  std::cout << "dump\n";

    for (std::size_t i = 0; i < cell.activ.size(); ++i) {
        outfile << '\n' << cell.activ[i].col.atomic;
        outfile << ' ' << cell.activ[i].vec[0];
        outfile << ' ' << cell.activ[i].vec[1];
        outfile << ' ' << cell.activ[i].vec[2];
    }

    for (std::size_t i = 0; i < cell.bound.size(); ++i) {
        outfile << '\n' << cell.bound[i].col.atomic + 99;
        outfile << ' ' << cell.bound[i].vec[0];
        outfile << ' ' << cell.bound[i].vec[1];
        outfile << ' ' << cell.bound[i].vec[2];
    }
    outfile << "\n";
}

void dump_supercell(Supercell const &cell) {
    static int frame = 0;

    static const std::string pre = "dump_";
    static const std::string post = ".xyz";

    dump_supercell(cell, pre + std::to_string(frame++) + post);
}

namespace {

// Parse "config["supercell"]["element_map"])", expects array of ["string", "species_string", "A/B"]
std::map<std::string, Colour> species_map_to_colour(
    toml::v2::table const &config,
    std::unordered_map<std::string, std::uint16_t> const &map) {
    // Key represents a name used to refer to atoms, value the corresponding colour

    if (toml::node_view types = config["supercell"]["element_map"]) {
        std::map<std::string, Colour> elems;
        //
        ALWAYS_CHECK(types.is_homogeneous(), "Element map incorrectly formatted");

        int i = 0;
        while (types[i]) {
            std::optional name = types[i][0].value<std::string>();
            std::optional species = types[i][1].value<std::string>();
            std::optional state = types[i][2].value<std::string_view>();

            ALWAYS_CHECK(name && species && state, "Element map incorrectly formatted");
            ALWAYS_CHECK(map.count(*species) == 1, "Mismatched element in configuration file");

            if (*state == "A") {
                elems.try_emplace(*name, Colour{map.at(*species), Colour::activ});
            } else if (*state == "B") {
                elems.try_emplace(*name, Colour{map.at(*species), Colour::bound});
            } else {
                ALWAYS_CHECK(false, "element-map only supports A/B options");
            }

            i++;
        }

        return elems;
    } else {
        ALWAYS_CHECK(false, "Configuration file missing \"supercell.element_map\"");
    }
}

}  // namespace

std::pair<Supercell, std::string> load_supercell(
    toml::v2::table const &config,
    std::unordered_map<std::string, std::uint16_t> const &species_map) {
    //
    std::pair<Supercell, std::string> out{
        Simbox::load(config),
        fetch<std::string>(config, "supercell", "in_file"),
    };

    // Parse xyz file
    std::ifstream file(out.second);

    ALWAYS_CHECK(file.good(), "Could not open supercell file");

    std::size_t num_atoms = 0;

    safe_getline(file) >> num_atoms;

    safe_getline(file);  // skip comment line

    std::map<std::string, Colour> elem_list = species_map_to_colour(config, species_map);

    for (size_t i = 0; i < num_atoms; ++i) {
        Vec3<double> vec;
        std::string atomic;

        safe_getline(file) >> atomic >> vec[0] >> vec[1] >> vec[2];

        auto it = elem_list.find(atomic);

        if (it != elem_list.end() && it->second.state == Colour::activ) {
            out.first.activ.emplace_back(out.first.canonicle_image(vec), it->second);
        } else if (it != elem_list.end() && it->second.state == Colour::bound) {
            out.first.bound.emplace_back(out.first.canonicle_image(vec), it->second);
        } else {
            ALWAYS_CHECK(false, "Unable to interpret input .xyz file:" + atomic);
        }
    }

    return out;
}
