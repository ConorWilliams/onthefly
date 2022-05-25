#pragma once

#include <cstddef>
#include <memory>
#include <utility>

#include "libatom/asserts.hpp"
#include "libatom/minimise/LBFGS/core.hpp"
#include "libatom/neighbour/neighbour_list.hpp"
#include "libatom/potentials/potential.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief A minimiser that uses the LBFGS algorithm.
   */
  class LBFGS {
  public:
    /**
     * @brief Used to configure the minimiser.
     */
    struct Options {
      /** @brief Number of previous steps held in memory. */
      std::size_t n = 10;
      /** @brief Number of steps before exit with failiure. */
      std::size_t iter_max = 2000;
      /** @brief Force convergence criterion (eV/Angstroms). */
      floating f2norm = 1e-5;
      /**
       * @brief Used to determine the skin size.
       *
       * Determined such that the expected number of neighbours is skin_frac * (num_neighbours if no
       * skin used). Hence must be larger than one. If larger then neighbour lists are built less
       * often but there will be more more non-neighbours in neighbour lists.
       */
      floating skin_frac = 1.1;
      /** @brief Trust tolerance, set larger to reduce trust radius change. */
      floating proj_tol = 0;
      /** @brief Maximum trust radius e.g max steps size (Angstroms). */
      floating max_trust = 0.5;
      /** @brief Minimum trust radius e.g initial step size (Angstroms). */
      floating min_trust = 0.05;
      /** @brief Trust radius expansion rate. */
      floating grow_trust = 1.5;
      /** @brief Trust radius contraction rate. */
      floating shrink_trust = 0.5;
      /** @brief Print out debug info and dumps minimisation trace to "lbfgs_debug.xyz" */
      bool debug = false;
    };

    explicit LBFGS(Options const &opt) : m_opt{opt}, m_core{opt.n} {}

    /**
     * @brief Move the atoms in the SimCell to a local minimum of the potential.
     *
     * @return true if minimisation converged.
     * @return false if failed to converge.
     */
    bool minimise(SimCell &, Potential &, std::size_t num_threads);

  private:
    Options m_opt;
    CoreLBFGS m_core;
    std::optional<NeighbourList> m_nl;
  };

}  // namespace otf