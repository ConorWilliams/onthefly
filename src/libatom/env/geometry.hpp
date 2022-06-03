#pragma once

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include "Eigen/src/Core/MatrixBase.h"
#include "fmt/core.h"
#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/utils.hpp"

namespace otf::env {
  /**
   * @brief Compute the centre of mass of a set of atom.
   */
  template <typename... M> [[nodiscard]] Vec3<floating> com(AtomVector<M...> const &ref) {
    Vec3<floating> sum = Vec3<floating>::Zero();

    for (auto const &elem : ref) {
      sum += elem(Position{});
    }

    return sum / ref.size();
  }

  /**
   * @brief Compute the rmsd between the position of the atoms in x and y after applying the
   * matrix transformation M to each atom in x.
   *
   * The rmsd (root mean squared distance) is the sqrt of the sum of the squared distances between
   * atoms.
   */
  template <typename... M1, typename... M2, typename E>
  [[nodiscard]] floating rmsd(E const &M, AtomVector<M1...> const &x, AtomVector<M2...> const &y) {
    //
    ASSERT(x.size() == y.size(), "Sizes must match!");

    floating sum_sq = 0;

    for (std::size_t i = 0; i < x.size(); ++i) {
      sum_sq += norm_sq(y[i](Position{}) - (M * x[i](Position{}).matrix()).array());
    }

    return std::sqrt(sum_sq);
  }

  /**
   * @brief Compute the rmsd between the positions of the atoms in x and y.
   *
   * The rmsd (root mean squared distance) is the sqrt of the sum of the squared distances between
   * atoms.
   */
  template <typename... M1, typename... M2>
  [[nodiscard]] floating rmsd(AtomVector<M1...> const &x, AtomVector<M2...> const &y) {
    return rmsd(Mat3<floating>::Identity(), x, y);
  }

  /**
   * @brief Compute the orthogonal transformation matrix that best transforms "x" onto "y".
   *
   * This does not permute any of the atoms it just minimizes the RMSD between the two sets.
   *
   * Uses the "Kabsch algorithm" see: https://en.wikipedia.org/wiki/Kabsch_algorithm
   */
  template <typename... M1, typename... M2>
  [[nodiscard]] Mat3<floating> ortho_onto(AtomVector<M1...> const &x, AtomVector<M2...> const &y) {
    //
    Mat3<floating> H = Mat3<floating>::Zero();

    ASSERT(x.size() == y.size(), "Can't rotate different sizes");

    for (std::size_t i = 0; i < x.size(); i++) {
      H += x[i](Position{}).matrix() * y[i](Position{}).matrix().transpose();
    }

    Eigen::JacobiSVD<Mat3<floating>> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);

    // Don't do sign correction as it's ok to reflect the coordinate system.

    return svd.matrixV() * svd.matrixU().transpose();
  }

  namespace detail {

    /** @brief The maximum number of atoms that can lie in the same plane. */
    inline constexpr std::size_t MAX_COPLANAR_ATOMS = 10;

    /**
     * @brief Helper function, verifies addition of n^th atom matches all previous atoms.
     */
    template <typename... M1, typename... M2> bool within_tol_up_to(AtomVector<M1...> const &mut,
                                                                    AtomVector<M2...> const &ref,
                                                                    floating tol, std::size_t n) {
      //
      for (std::size_t i = 0; i < std::min(n, MAX_COPLANAR_ATOMS); ++i) {
        // Intra atomic distances
        floating ref_ni = norm(ref[n](Position{}) - ref[i](Position{}));
        floating mut_ni = norm(mut[n](Position{}) - mut[i](Position{}));

        if (std::abs(ref_ni - mut_ni) > tol) {
          return false;
        }
      }
      return true;
    }

    /**
     * @brief Recursive implementation of for_equiv_perms.
     *
     * Uses the GREEDY_PERM algorithm, recursively finds a new atom, such that the order that
     * the intra-atom distances in each AtomVector match (within tolerance SQRT_2 * delta). The
     * final permutation must be able to be transformed (via ortho_onto(mut, other)) such that the
     * RMSD between the two sets is less than delta.
     */
    template <typename... M1, typename... M2, typename F>
    bool for_equiv_perms_impl(AtomVector<M1...> &mut, AtomVector<M2...> const &ref, floating delta,
                              std::size_t n, F &&f) {
      // Termination criterion.
      if (n >= mut.size()) {
        //
        Mat3<floating> O = ortho_onto(mut, ref);

        if (floating dr = rmsd(O, mut, ref); dr < delta) {
          return std::invoke(f, O, dr);
        } else {
          return false;
        }
      }

      using std::swap;  // ADL

      // Attempt to find an atom to put in i^th position.
      for (std::size_t i = n; i < ref.size(); ++i) {
        if (mut[i](Colour{}) == ref[n](Colour{})) {
          // Test next candidate
          swap(mut[n], mut[i]);

          // Verify all other distances then recurse
          if (within_tol_up_to(mut, ref, delta * M_SQRT2, n)) {
            if (for_equiv_perms_impl(mut, ref, delta, n + 1, f)) {
              // Enable early exit.
              return true;
            }
          }

          // Undo swap such that if we have to back-track order is the same.
          swap(mut[n], mut[i]);
        }
      }

      return false;
    }

  }  // namespace detail

  /**
   * @brief Create a function object to explore equivalent permutations of another AtomVector.
   *
   * This function is for greedily exploring the permutations of atoms in another AtomVector such
   * that, after an orthogonal transformation generated by ortho_onto(), the rmsd between the two
   * sets of atoms is less than delta.
   *
   * Example:
   *
   * @code{.cpp}
   *
   * for_equiv_perms(ref, 0.2)(mut, [](Mat3<floating> const & O, floating rmsd){
   *     // If this function is called then "mut" has been permuted into an equivalent permutation.
   *     // "O" is the matrix that maps "mut" to "ref" that "rmsd" == rmsd(O, mut, ref).
   *
   *     // We could now do something with "rmsd" and "O". We must not mutate "ref" or "mut".
   *
   *     // If we "return true" then exploration/function terminates.
   *     // If we "return false" then exploration/function continues.
   * });
   *
   * @endcode
   *
   */
  template <typename... M> auto for_equiv_perms(AtomVector<M...> const &ref, floating delta) {
    return [&ref, delta](auto &mut, auto &&callback) {
      ASSERT(ref.size() == mut.size(), "Sizes must match!");

      if constexpr (std::is_same_v<decltype(mut), AtomVector<M...>>) {
        ASSERT(&ref != &mut, "Cannot perm onto self.");
      }

      detail::for_equiv_perms_impl(mut, ref, delta, 0, callback);
    };
  }

  /**
   * @brief A geometry models a local distribution of periodically resolved Atoms.
   *
   * The distribution is centred on the first atom.
   */
  template <typename... Mems> class Geometry : public AtomVector<Mems...> {
  public:
    /**
     * @brief Returned by permute_onto(delta, x), contains extra info about the permutation.
     */
    struct PermResult {
      /** @brief  Orthogonal transformation required to map this onto x. */
      Mat3<floating> O;
      /** @brief Equal to rmsd(O, *this, x). */
      floating rmsd;
    };

    /**
     * @brief Attempt to permute the atoms in this Geometry into the same order as the atoms in
     * other.
     *
     * @return std::optional<PermResult> An engaged PermResult if a permutation was found.
     */
    template <typename... Ts>
    [[nodiscard]] std::optional<PermResult> permute_onto(AtomVector<Ts...> const &other,
                                                         floating delta) {
      std::optional<PermResult> res = std::nullopt;

      for_equiv_perms(other, delta)(*this, [&](Mat3<floating> const &O, floating rmsd) {
        res = PermResult{O, rmsd};
        // First match accepted
        return true;
      });

      return res;
    }

  private:
  };

}  // namespace otf::env