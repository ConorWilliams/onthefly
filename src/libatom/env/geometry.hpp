#pragma once

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <cmath>
#include <cstddef>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include "fmt/core.h"
#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/utils.hpp"

namespace otf::env {

  /** @brief The maximum number of atoms that can lie in the same plane. */
  inline constexpr std::size_t MAX_COPLANAR_ATOMS = 10;

  /**
   * @brief A geometry models a local distribution of periodically resolved Atoms.
   *
   * The distribution is centred on the first atom.
   */
  template <typename... Mems> class Geometry : public AtomVector<Mems...> {
  public:
    /**
     * @brief Compute the centre of mass of the atoms in the geometry
     */
    [[nodiscard]] Vec3<floating> com() const {
      Vec3<floating> sum = Vec3<floating>::Zero();

      for (auto const &elem : *this) {
        sum += elem(Position{});
      }

      return sum / this->size();
    }

    /**
     * @brief Compute the root mean squared distance between the position of the atoms in this and
     * x.
     *
     * This is the sqrt of the sum of the squared distances between atoms.
     */
    template <typename... Ts> [[nodiscard]] floating rmsd(AtomVector<Ts...> const &x) const {
      floating sum = 0;
      for (std::size_t i = 0; i < this->size(); i++) {
        sum += norm_sq((*this)[i](Position{}) - x[i](Position{}));
      }
      return std::sqrt(sum);
    }

    /**
     * @brief Compute the orthogonal transformation matrix that best transforms this onto "other".
     *
     * This does not permute any of the atoms it just minimizes the RMSD between the two sets.
     *
     * Uses the "Kabsch algorithm" see: https://en.wikipedia.org/wiki/Kabsch_algorithm
     */
    template <typename... Ts>
    [[nodiscard]] Mat3<floating> ortho_onto(AtomVector<Ts...> const &other) const {
      //
      Mat3<floating> H = Mat3<floating>::Zero();

      ASSERT(this->size() == other.size(), "Can't rotate different sizes");

      for (std::size_t i = 0; i < this->size(); i++) {
        H += (*this)[i](Position{}).matrix() * other[i](Position{}).matrix().transpose();
      }

      Eigen::JacobiSVD<Mat3<floating>> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);

      // Don't do sign correction as it's ok to reflect coordinate system

      return svd.matrixV() * svd.matrixU().transpose();
    }

    /**
     * @brief Returned by .permute_onto(delta, x) that contains extra info about the permutation.
     */
    struct PermResult {
      /** @brief Minimum RMSD between this and x */
      floating rmsd;
      /** @brief  Orthogonal transformation required to map this onto x. */
      Mat3<floating> O;
    };

    /**
     * @brief Attempt to permute the atoms in this Geometry into the same order as the atoms in
     * other.
     *
     * Uses the GREEDY_PERM algorithm, matches the order such that the intra-atom distances in each
     * set match. The final permutation must be able to be transformed (via .ortho_onto(other)) such
     * that the RMSD between the two sets is less than delta.
     *
     * @return std::optional<PermResult> An engaged PermResult if a permutation was found.
     */
    template <typename... Ts>
    [[nodiscard]] std::optional<PermResult> permute_onto(AtomVector<Ts...> const &other,
                                                         floating delta) {
      //
      ASSERT(this->size() == other.size(), "Sizes must match!");

      if constexpr (std::conjunction_v<std::is_same<Ts, Mems>...>) {
        ASSERT(static_cast<AtomVector<Ts...> *>(this) != &other, "Cannot perm onto self.");
      }

      return _permute_onto(*this, other, delta, 1);
    }

  private:
    /**
     * @brief Recursive implementation of perm_onto.
     */
    template <typename... Ts>
    static std::optional<PermResult> _permute_onto(Geometry<Mems...> &mut,
                                                   AtomVector<Ts...> const &ref, floating delta,
                                                   std::size_t n) {
      // Termination criterion
      if (n >= mut.size()) {
        Mat3<floating> R = mut.ortho_onto(ref);

        floating sum_sq = 0;

        for (std::size_t i = 0; i < mut.size(); ++i) {
          sum_sq += norm_sq(ref[i](Position{}) - (R * mut[i](Position{}).matrix()).array());
        }

        if (sum_sq < delta * delta) {
          fmt::print("Det={}, Trace={} ", R.determinant(), R.trace());

          if (R.determinant() < 0) {
            fmt::print("Theta={:<3}deg : ", 360. / 2 / M_PI * std::acos(0.5 * (R.trace() + 1)));
          } else {
            fmt::print("Theta={:<3}deg : ", 360. / 2 / M_PI * std::acos(0.5 * (R.trace() - 1)));
          }

          for (std::size_t i = 0; i < mut.size(); i++) {
            fmt::print("{}, ", mut[i](Index{}));
          }

          fmt::print("\n");

          return PermResult{std::sqrt(sum_sq), R};
        } else {
          return std::nullopt;
        }
      }

      using std::swap;  // ADL

      // Attempt to find an atom to put in i^th position
      for (std::size_t i = n; i < ref.size(); ++i) {
        if (mut[i](Colour{}) == ref[n](Colour{})) {
          // Test next candidate
          swap(mut[n], mut[i]);

          // Verify all other distances then recurse
          if (within_tol_up_to(mut, ref, delta * M_SQRT2, n)) {
            if (std::optional rot = _permute_onto(mut, ref, delta, n + 1)) {
              // return rot;
            }
          }

          // Undo swap such that if we have to back-track order is the same
          swap(mut[n], mut[i]);
        }
      }

      return std::nullopt;
    }

    /**
     * @brief Helper function for permute_onto, verifies addition of n^th atom matches all previous
     * atoms.
     */
    template <typename... Ts> static bool within_tol_up_to(Geometry<Mems...> &mut,
                                                           AtomVector<Ts...> const &ref,
                                                           floating tol, std::size_t n) {
      //  for (std::size_t i = 0; i < n; ++i)
      //  for (std::size_t i = 0; i < std::min(n, MAX_COPLANAR_ATOMS); ++i)
      for (std::size_t i = 0; i < std::min(n, MAX_COPLANAR_ATOMS); ++i) {
        if (std::abs(norm(ref[n](Position{}) - ref[i](Position{}))
                     - norm(mut[n](Position{}) - mut[i](Position{})))
            > tol) {
          return false;
        }
      }

      return true;
    }
  };

}  // namespace otf::env