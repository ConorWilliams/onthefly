
#include "package/package.hpp"

#include <stdexcept>

#include "config.hpp"
#include "local/environment.hpp"
#include "supercell.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

namespace options {

Packager Packager::load(toml::v2::table const &config) {
    Packager opt;

    opt.mode = config["package"]["mode"].value_or(opt.mode);

    ALWAYS_CHECK(std::find(opt.valid.begin(), opt.valid.end(), opt.mode) != opt.valid.end(),
                 "Invalid packaging mode");

    opt.unpack_tol = fetch<double>(config, "package", "unpack_tol");

    if (opt.mode != "global") {
        opt.r_active = fetch<double>(config, "package", "r_active");
        opt.r_boundary = fetch<double>(config, "package", "r_boundary");
        opt.require_centre = fetch<bool>(config, "package", "require_centre");
    }

    opt.mech = Mechanism::load(config);

    return opt;
}

}  // namespace options

Package Packager::pack_full(Supercell const &cell, std::size_t centre) {
    // Package to return
    Package pkg(cell.activ.size());

    static_cast<Supercell &>(pkg.subcell) = cell;

    pkg.subcell.centre = centre;

    return pkg;
}

Package Packager::pack_unresolved(Supercell const &cell, std::size_t centre) {
    // Package to return
    Package pkg(cell.activ.size());

    // Unresolved subcell has same simbox as initial supercell
    static_cast<Simbox &>(pkg.subcell) = static_cast<Simbox const &>(cell);

    // TODO : could make fwd_map a member of catalogue and reuse

    // Mapping from Subcell -> Supercell
    std::vector<std::size_t> fwd_map;

    Vec3<double> cen = cell.activ[centre].vec;

    // Populate subcell from activ;
    for (std::size_t i = 0; i < cell.activ.size(); i++) {
        //
        double dr = norm_sq(cell.min_image(cell.activ[i].vec, cen));

        if (dr < _opt.r_active * _opt.r_active) {
            fwd_map.push_back(i);
            pkg.subcell.activ.emplace_back(cell.activ[i].vec, cell.activ[i].col);
        } else if (dr < _opt.r_boundary * _opt.r_boundary) {
            pkg.subcell.bound.emplace_back(cell.activ[i].vec, cell.activ[i].col);
        }
    }

    // Populate subcell from bound;
    for (size_t i = 0; i < cell.bound.size(); i++) {
        double dr = norm_sq(cell.min_image(cen, cell.bound[i].vec));

        if (dr < _opt.r_boundary * _opt.r_boundary) {
            pkg.subcell.bound.emplace_back(cell.bound[i].vec, cell.bound[i].col);
        }
    }

    // Build reverse map: Supercell -> Subcell
    for (std::size_t i = 0; i < fwd_map.size(); i++) {
        pkg.map[fwd_map[i]] = i;
    }

    ALWAYS_CHECK(pkg.map[centre], "Centre not included in subcell");

    pkg.subcell.centre = *pkg.map[centre];

    return pkg;
}

std::vector<Package> Packager::pack(Supercell const &cell, std::vector<std::size_t> centre) {
    //
    std::vector<Package> out;

    if (_opt.mode == "global") {
        for (auto &&cen : centre) {
            out.push_back(pack_full(cell, cen));
        }
    } else if (_opt.mode == "local") {
        for (auto &&cen : centre) {
            out.push_back(pack_unresolved(cell, cen));
        }
    } else {
        throw std::runtime_error("Packaging mode \"" + _opt.mode + "\" invalid");
    }

    return out;
}

/////////////////////////////////////////////////////////////////////////////////////////

namespace {

// Returns the absolute and fraction of the mechanisms captured during localisation.
std::pair<double, double> quality(ProtoMech const &proto, std::vector<Vec3<double>> const &disp) {
    double cap = [&] {
        double sum = 0;

        for (auto const &elem : disp) {
            sum += norm_sq(elem);
        }
        return std::sqrt(sum);
    }();

    double tot = norm(proto.disp);

    return {cap, cap / tot};
}

// Convert proto-mechanisms to local mechanism, locality determined by reference geometry, rotor and
// mapping, only active atoms in the geometry have their motion recorded.
std::optional<Mechanism> localise(ProtoMech const &proto,
                                  Geometry const &ref,
                                  Mat3<double> const &rotor,
                                  std::vector<std::optional<std::size_t>> const &map) {
    std::vector<Vec3<double>> disp;

    // Rearrange mechanism onto canonical order
    for (std::size_t i = 0; i < ref.size(); i++) {
        if (ref[i].col.state == Colour::activ) {
            if (map[ref[i].idx]) {
                // Atom index in subcell
                std::size_t j = *map[ref[i].idx];

                Vec3<double> delta = {
                    proto.disp[3 * j + 0],
                    proto.disp[3 * j + 1],
                    proto.disp[3 * j + 2],
                };

                disp.emplace_back(rotor * delta.matrix());
            } else {
                return std::nullopt;
            }
        }
    }

    auto [abs, rel] = quality(proto, disp);

    return Mechanism{proto, std::move(disp), abs, rel};
}

double l2_similarity(Geometry const &geo, Mat3<double> const &R, Geometry const &ref) {
    double sum_sq = 0;

    for (std::size_t i = 0; i < geo.size(); ++i) {
        sum_sq += norm_sq(ref[i].vec - (R * geo[i].vec.matrix()).array());
    }

    return std::sqrt(sum_sq);
}

}  // namespace

void Packager::unpack_full(Package &&pkg,
                           std::vector<Geometry> const &geos,
                           std::vector<Catalogue::pointer> const &env) {
    for (auto &&proto : pkg.mechs) {
        // Centre in frame of subcell
        std::size_t centre = proto.find_centre();

        Mat3<double> R = geos[centre].rotor_onto(env[centre]->geo);

        // Verify centre and env are exceptional match
        if (l2_similarity(geos[centre], R, env[centre]->geo) > _opt.unpack_tol) {
            continue;
        }

        std::vector<Vec3<double>> disp;

        // Rearrange mechanism onto canonical order
        for (std::size_t i = 0; i < geos[centre].size(); i++) {
            // Ignore inactive elements of geometry
            if (geos[centre][i].col.state == Colour::activ) {
                Vec3<double> delta = {
                    proto.disp[3 * geos[centre][i].idx + 0],
                    proto.disp[3 * geos[centre][i].idx + 1],
                    proto.disp[3 * geos[centre][i].idx + 2],
                };

                disp.emplace_back(R * delta.matrix());
            }
        }

        auto [abs, rel] = quality(proto, disp);

        env[centre]->try_push_mech({proto, std::move(disp), abs, rel},  // Construct mech
                                   _opt.mech.energy_abs_tol,
                                   _opt.mech.energy_frac_tol,
                                   _opt.mech.r_tol);
    }
}

void Packager::unpack_unresolved(Package &&pkg,
                                 std::vector<Geometry> const &geos,
                                 std::vector<Catalogue::pointer> const &env) {
    for (auto &&proto : pkg.mechs) {
        // Centre in frame of subcell
        std::size_t mech_centre = proto.find_centre();

        if (!_opt.require_centre || mech_centre == pkg.subcell.centre) {
            // Centre in frame of supercell
            std::size_t super_centre = [&] {
                for (std::size_t i = 0; i < pkg.map.size(); i++) {
                    if (pkg.map[i] == mech_centre) {
                        return i;
                    }
                }
                ALWAYS_CHECK(false, "Could not find centre in pkg.map");
            }();

            Mat3<double> R = geos[super_centre].rotor_onto(env[super_centre]->geo);

            if (l2_similarity(geos[super_centre], R, env[super_centre]->geo) > _opt.unpack_tol) {
                continue;
            }

            if (std::optional mech = localise(proto, geos[super_centre], R, pkg.map)) {
                env[super_centre]->try_push_mech(std::move(*mech),
                                                 _opt.mech.energy_abs_tol,
                                                 _opt.mech.energy_frac_tol,
                                                 _opt.mech.r_tol);
            } else {
                std::cerr << "WARNING: subcell does not envelope geometry, .require_centre="
                          << _opt.require_centre << "\n";
            }
        } else {
            // std::cout << "Mech centre does not align\n";
        }
    }
}

void Packager::unpack(std::vector<Package> &&pkgs,
                      std::vector<Geometry> const &geos,
                      std::vector<Catalogue::pointer> const &env) {
    if (_opt.mode == "global") {
        for (auto &&pkg : pkgs) {
            unpack_full(std::move(pkg), geos, env);
        }
    } else if (_opt.mode == "local") {
        for (auto &&pkg : pkgs) {
            unpack_unresolved(std::move(pkg), geos, env);
        }
    } else {
        throw std::runtime_error("Packaging mode \"" + _opt.mode + "\" invalid");
    }
}