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
     * @brief Get a const reference to the discrete key.
     */
    Key const &key() const noexcept { return m_key; }

    /**
     * @brief Rebuild this local environment to be the env of the idx'th atom.
     */
    void rebuild(SimCell const &atoms, neighbour::List const &nl, std::size_t idx);

  private:
    /**
     * @brief An ordered representation of the inta-atomic distances in the Topology
     */
    class Fingerprint {
    public:
      /**
       * @brief A fast test to see if two local environemnts **may** be equivilant.
       *
       * @param a
       * @param b
       * @param delta
       * @return true
       * @return false
       */
      friend bool equiv(Fingerprint const &a, Fingerprint const &b, double delta) {
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

      void clear() noexcept {
        m_r_0j.clear();
        m_r_ij.clear();
      }

    private:
      friend class LocalEnv;

      std::vector<floating> m_r_0j;
      std::vector<floating> m_r_ij;
    };

    Key m_key;

    Fingerprint m_fingerprint;
  };

}  // namespace otf::env