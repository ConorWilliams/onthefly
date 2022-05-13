#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "cereal/types/vector.hpp"
#include "config.hpp"
#include "fuzzy_key.hpp"
#include "supercell.hpp"
#include "utility.hpp"

// Class representing the real space, periodicity resolved, representation of a (sub)set of
// atoms in a supercell. Embeds the fuzzy_key associated with the atoms in the geometry.
class Geometry {
  private:
    struct geo_atom {
        Vec3<double> vec;
        Colour col;
        std::size_t idx;

        geo_atom(Vec3<double> const &vec, Colour const &col, std::size_t idx)
            : vec(vec), col(col), idx(idx) {}

        geo_atom() = default;

        // Deliberate non-serialisation of .idx
        template <class Archive> void serialize(Archive &ar) { ar(vec, col); }
    };

  public:
    // Expose minimal vector-like API
    std::size_t size() const { return _atoms.size(); }

    geo_atom &operator[](std::size_t i) { return _atoms[i]; }
    geo_atom const &operator[](std::size_t i) const { return _atoms[i]; }

    void clear() {
        _atoms.clear();
        _fuzzy_key.clear();
    }

    template <typename... Args> void emplace_back(Args &&...args) {
        _atoms.emplace_back(std::forward<Args>(args)...);
    }

    // Must be called after all atoms emplaced
    void finalise();

    // Return the transformation matrix R that transforms this onto "other". Uses the "Kabsch
    // algorithm" computes optimum rotation see: https://en.wikipedia.org/wiki/Kabsch_algorithm
    Mat3<double> rotor_onto(Geometry const &other) const;

    struct Result {
        double dr;       // l2-norm
        Mat3<double> R;  // Required rotation
    };

    // Attempts to permute atoms in this such that after rotation by the returned matrix R the l2
    // norm of the distance between the two point sets is dr < delta
    std::optional<Result> permute_onto(double delta, Geometry const &ref);

    double norm(Geometry const &other) const;

    bool equiv(double tol, Geometry const &other) const;

    template <class Archive> void serialize(Archive &ar) { ar(_atoms, _fuzzy_key); }

  private:
    std::vector<geo_atom> _atoms{};
    fuzzy_key _fuzzy_key;
};
