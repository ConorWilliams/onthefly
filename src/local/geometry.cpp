#include "geometry.hpp"

#include <cmath>
#include <optional>
#include <utility>
#include <vector>

#include "cereal/access.hpp"
#include "cereal/types/vector.hpp"
#include "config.hpp"
#include "local/fuzzy_key.hpp"
#include "utility.hpp"

// Must be called after all atoms emplaced
void Geometry::finalise() {
    // Sanity checks
    CHECK(size() > 0, "Too few atoms in geometry");

    // Make COM {0, 0, 0}
    Vec3<double> sum = {0, 0, 0};

    for (auto const &elem : _atoms) {
        sum += elem.vec;
    }

    Vec3<double> const com = sum / _atoms.size();

    for (auto &&elem : _atoms) {
        elem.vec -= com;
    }

    // Order the atoms in a way that accelerates calls to permute_onto(), first atom (centre) is not
    // moved.
    std::sort(std::next(_atoms.begin()), _atoms.end(), [&](auto const &a, auto const &b) {
        if (a.col != b.col) {
            return a.col < b.col;
        } else {
            return norm_sq(a.vec) < norm_sq(b.vec);
        }
    });

    // Build _fuzzy_key
    _fuzzy_key.build(_atoms);
}

// Return the transformation matrix R that transforms this onto "other". Uses the "Kabsch
// algorithm" computes optimum rotation see: https://en.wikipedia.org/wiki/Kabsch_algorithm
Mat3<double> Geometry::rotor_onto(Geometry const &other) const {
    Mat3<double> H = Mat3<double>::Zero();

    CHECK(size() == other.size(), "Can't rotate different sizes");

    for (std::size_t k = 0; k < size(); ++k) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                H(i, j) += _atoms[k].vec[i] * other[k].vec[j];
            }
        }
    }

    Eigen::JacobiSVD<Mat3<double>> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);

    // Don't do sign correction as it's ok to reflect coordinate system

    return svd.matrixV() * svd.matrixU().transpose();
}

double Geometry::norm(Geometry const &other) const {
    return fuzzy_norm(_fuzzy_key, other._fuzzy_key);
}

bool Geometry::equiv(double tol, Geometry const &other) const {
    return equivalent(tol, _fuzzy_key, other._fuzzy_key);
}

namespace {

// Helper function for permute_onto, verifies addition of n^th atom matches all previous atoms
bool within_tol_up_to(Geometry const &ref, Geometry const &mut, double tol, std::size_t n) {
    //  for (std::size_t i = 0; i < n; ++i)
    //  for (std::size_t i = 0; i < std::min(n, MAX_COPLANAR_ATOMS); ++i)
    for (std::size_t i = 0; i < std::min(n, MAX_COPLANAR_ATOMS); ++i) {
        if (std::abs(norm(ref[n].vec - ref[i].vec) - norm(mut[n].vec - mut[i].vec)) > tol) {
            return false;
        }
    }

    return true;
}

std::optional<Geometry::Result> permute_onto(Geometry const &ref,
                                             Geometry &mut,
                                             double delta,
                                             std::size_t n) {
    // Termination criterion
    if (n >= mut.size()) {
        Mat3<double> R = mut.rotor_onto(ref);

        double sum_sq = 0;

        for (std::size_t i = 0; i < mut.size(); ++i) {
            sum_sq += norm_sq(ref[i].vec - (R * mut[i].vec.matrix()).array());
        }

        if (sum_sq < delta * delta) {
            return Geometry::Result{std::sqrt(sum_sq), R};
        } else {
            return std::nullopt;
        }
    }

    using std::swap;  // ADL

    // Attempt to find an atom to put in i^th position
    for (std::size_t i = n; i < ref.size(); ++i) {
        if (mut[i].col == ref[n].col) {
            // Test next candidate
            swap(mut[n], mut[i]);

            // Verify all other distances then recurse
            if (within_tol_up_to(ref, mut, delta * SQRT_2, n)) {
                if (auto rot = permute_onto(ref, mut, delta, n + 1)) {
                    return rot;
                }
            }

            // Undo swap such that if we have to back-track order is the same
            swap(mut[n], mut[i]);
        }
    }

    return std::nullopt;
}

}  // namespace

std::optional<Geometry::Result> Geometry::permute_onto(double delta, Geometry const &ref) {
    CHECK(this != &ref, "Cannot perm onto self");
    CHECK(size() == ref.size(), "Wrong number of atoms in other");
    return ::permute_onto(ref, *this, delta, 1);
}
