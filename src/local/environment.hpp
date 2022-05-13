#pragma once

#include <optional>
#include <vector>

#include "cereal/types/base_class.hpp"
#include "cereal/types/vector.hpp"
#include "config.hpp"
#include "local/geometry.hpp"
#include "supercell.hpp"
#include "utility.hpp"

namespace options {

struct Mechanism {
    double r_tol;            // (Angstroms) L2 between mechanisms to be considered distinct
    double energy_abs_tol;   // (eV) For mechanisms to be considered distinct
    double energy_frac_tol;  // For mechanisms to be considered distinct
    double rel_cap_tol;      // Fraction of mechanism required to be captured by LE

    static Mechanism load(toml::v2::table const& config);
};

}  // namespace options

class MechBase {
  public:
    double activ_energy;
    double delta_energy;
    double pre_factor;

    template <class Archive> void serialize(Archive& ar) {
        ar(activ_energy, delta_energy, pre_factor);
    }

  protected:
    bool within_tol(MechBase const& other, double abs_tol, double frac_tol) const;
};

// Stores the displacement vectors for each atom in a geometry during some mechanism alongside rate
// metadata
class Mechanism : public MechBase {
  public:
    std::vector<Vec3<double>> disp;

    double abs_cap;  // L2 norm of the mechanisms displacement
    double rel_cap;  // Fraction of the L2 norm of the ProtoMech captured by this mechanism

    bool within_tol(Mechanism const& other, double abs_tol, double frac_tol, double r_tol) const;

    template <class Archive> void serialize(Archive& ar) {
        ar(cereal::base_class<MechBase>(this), disp, abs_cap, rel_cap);
    }
};

// Strongly typed mechanism that has not been refined to a local environment.
class ProtoMech : public MechBase {
  public:
    VecN<double> disp;

    bool within_tol(ProtoMech const& other, double abs_tol, double frac_tol, double r_tol) const;

    std::size_t find_centre() const;
};

// Represents a topological basin, stores mechanisms escaping basin
struct Environment {
    Geometry geo;                    // Reference geometry (pre_ordered)
    double delta;                    //
    int freq = 0;                    // Occurrence count
    int search = 0;                  // Number of sps initiated from topology
    std::vector<Mechanism> mechs{};  // All mechanisms from this basin

    Environment(Geometry const& geo, double del) : geo{geo}, delta(del) {}

    Environment() = default;  // Cereal

    bool try_push_mech(Mechanism&& m, double abs_tol, double frac_tol, double r_tol);

    template <class Archive> void serialize(Archive& ar) { ar(geo, delta, freq, search, mechs); }
};
