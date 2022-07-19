#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/env/geometry.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/utils.hpp"

namespace otf::env {

  class LocalEnv;

  /**
   * @brief A local environment represent a localised geometry augemented with a key and
   * fingerprint.
   */
  class LocalEnv : public Geometry<Position, Colour, Index> {
  public:
    /**
     * @brief The discrete representation of this environment.
     *
     * The first pair is the colour of the centre atom (for .first and .second = 1), the remaining
     * pairs are .first = colour and .second = count.
     */
    using Key = std::vector<std::pair<std::size_t, std::size_t>>;

    /**
     * @brief An ordered representation of the inta-atomic distances in this Geometry.
     */
    class Fingerprint {
    public:
      /**
       * @brief A fast test to see if two local environemnts **may** be equivilant.
       *
       * Explicitly this tests that the paired ordered intra-atomic distances are within "delta".
       * If using this as a filter for permute_** algorithims then using "delta" = M_SQRT2 *
       * delta_for_perm will
       */
      friend bool equiv(Fingerprint const& a, Fingerprint const& b, double delta) {
        if (a.m_r_0j.size() != b.m_r_0j.size()) {
          return false;
        }

        for (std::size_t i = 0; i < a.m_r_0j.size(); i++) {
          if (std::abs(a.m_r_0j[i] - b.m_r_0j[i]) > delta) {
            return false;
          }
        }

        ASSERT(a.m_r_ij.size() == b.m_r_ij.size(), "Secondary distances should match");

        for (std::size_t i = 0; i < a.m_r_ij.size(); i++) {
          if (std::abs(a.m_r_ij[i] - b.m_r_ij[i]) > delta) {
            return false;
          }
        }

        return true;
      }

      floating r_min() const {
        ASSERT(!m_r_0j.empty() && !m_r_ij.empty(), "Not enough atoms for f_min!");
        return std::min(m_r_0j[0], m_r_ij[0]);
      }

    private:
      friend class LocalEnv;

      void clear() noexcept {
        m_r_0j.clear();
        m_r_ij.clear();
      }

      std::vector<floating> m_r_0j;
      std::vector<floating> m_r_ij;
    };

    /**
     * @brief Get a const reference to the discrete key.
     */
    Key const& key() const noexcept { return m_key; }

    /**
     * @brief Get a const reference to the fingerprint.
     */
    Fingerprint const& fingerprint() const noexcept { return m_fingerprint; }

    /**
     * @brief Rebuild this local environment to be the env of the idx'th atom.
     */
    void rebuild(SimCell const& atoms, neighbour::List const& nl, std::size_t idx);

  private:
    Key m_key;

    Fingerprint m_fingerprint;
  };

  /**
   * @brief This class builds and stores a list of local environments from a SimCell.
   */
  class EnvCell {
  public:
    struct Options {
      /** @brief Radius of environment. */
      floating r_env;
    };

    /**
     * @brief Construct a new Env Cell object.
     */
    EnvCell(Options const& opt, OrthoSimBox const& box) : m_nl(box, opt.r_env){};

    /**
     * @brief Get the ith local environment
     */
    LocalEnv const& env(std::size_t i) const noexcept { return m_envs[i]; }

    /**
     * @brief Get the ith local environment
     */
    LocalEnv& env(std::size_t i) noexcept { return m_envs[i]; }

    /**
     * @brief Construct the local environemnts around all the atom.
     */
    void rebuild(SimCell const& atoms, std::size_t num_threads);

  private:
    std::vector<LocalEnv> m_envs;

    neighbour::List m_nl;
  };

}  // namespace otf::env