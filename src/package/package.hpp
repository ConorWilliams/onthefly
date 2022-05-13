#pragma once

#include <future>

#include "config.hpp"
#include "local/catalogue.hpp"
#include "local/environment.hpp"
#include "supercell.hpp"
#include "utility.hpp"

namespace options {

struct Packager {
    static constexpr std::array valid = {"global", "local"};

    std::string mode = "global";  // Must be in valid

    double unpack_tol;    // (Angstroms) maximum l2 deviation between env and ref to unpack
    double r_active;      // (Angstroms), Radius in which atoms can move during SPS
    double r_boundary;    // (Angstroms), Radius for inclusion in subcell as bound atoms
    bool require_centre;  // Require mechanisms centred on displacement centre

    Mechanism mech;

    static Packager load(toml::v2::table const &config);
};

}  // namespace options

class Packager;

// Packages work, results space and metadata required to process results
class Package {
  public:
    Workcell subcell;                             // Cell in which to conduct SPS
    std::future<std::vector<ProtoMech>> f_mechs;  // Future to store resulting mechanisms
    std::vector<ProtoMech> mechs;

  private:
    friend class Packager;

    explicit Package(std::size_t active_size) : map(active_size, std::nullopt) {}

    std::vector<std::optional<std::size_t>> map;  // Perform mapping from [super->sub]cell idxs
};

class Packager {
  public:
    explicit Packager(options::Packager const &opt) : _opt(opt) {}

    // Make a Package centred on each atom in centres
    std::vector<Package> pack(Supercell const &cell, std::vector<std::size_t> centres);

    // Unpack a set of proto-mechanisms (in a package) into the corresponding local environments
    void unpack(std::vector<Package> &&pkgs,
                std::vector<Geometry> const &geos,
                std::vector<Catalogue::pointer> const &env);

  private:
    options::Packager _opt;

    Package pack_full(Supercell const &cell, std::size_t centre);
    Package pack_unresolved(Supercell const &cell, std::size_t centre);

    void unpack_full(Package &&pkgs,
                     std::vector<Geometry> const &geos,
                     std::vector<Catalogue::pointer> const &env);

    void unpack_unresolved(Package &&pkgs,
                           std::vector<Geometry> const &geos,
                           std::vector<Catalogue::pointer> const &env);
};
