#pragma once

#include <array>
#include <cstddef>
#include <fstream>

#include "libatom/asserts.hpp"
#include "libatom/potentials/spline.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief Holds the parsed tabulated EAM data and reconstruct f,phi,v smoothly via splines.
   */
  class DataEAM {
  public:
    /**
     * @brief Parse the lammps eam/fs style stream and construct an DataEAM object.
     *
     * For detail on file specification see: https://docs.lammps.org/pair_eam.html
     */
    DataEAM(std::ifstream in);

    floating rcut() const { return m_rcut; }

    /**
     * @brief Fetch the embedding energy function corresponding to the atom with atomic number 'a'.
     */
    Spline const &f(std::size_t a) const {
      ASSERT(a < 112, "Invalid atomic number");
      ASSERT(m_atomic2idx[a] < m_num_species, "Species not in this potential");
      return m_f[m_atomic2idx[a]];
    }

    /**
     * @brief Fetch the electron density function corresponding to the atom pair with atomic numbers
     * 'a' and 'b'.
     */
    Spline const &phi(std::size_t a, std::size_t b) const {
      ASSERT(a < 112 && b < 112, "Invalid atomic number");
      ASSERT(m_atomic2idx[a] < m_num_species, "Species not in this potential");
      ASSERT(m_atomic2idx[b] < m_num_species, "Species not in this potential");
      return m_phi[index(m_atomic2idx[a], m_atomic2idx[b])];
    }

    /**
     * @brief Fetch the symetric pair potential function corresponding to the atom pair with atomic
     * numbers 'a' and 'b'.
     */
    Spline const &v(std::size_t a, std::size_t b) const {
      ASSERT(a < 112 && b < 112, "Species out of bounds");
      ASSERT(m_atomic2idx[a] < m_num_species, "Species not in this potential");
      ASSERT(m_atomic2idx[b] < m_num_species, "Species not in this potential");
      return m_v[sym_index(m_atomic2idx[a], m_atomic2idx[b])];
    }

    /**
     * @brief Fetch the mass of the atom with atomic number 'a'.
     */
    floating const &mass(std::size_t a) const {
      ASSERT(a < m_num_species, "Invalid atomic number");
      ASSERT(m_atomic2mass[a] != 0, "Species not in this potential");
      return m_atomic2mass[a];
    }

  private:
    std::size_t m_num_species;

    floating m_rcut;

    std::array<std::size_t, 112> m_atomic2idx;
    std::array<floating, 112> m_atomic2mass;

    std::vector<Spline> m_f;
    std::vector<Spline> m_phi;
    std::vector<Spline> m_v;

    std::size_t index(std::size_t i, std::size_t j) const {
      ASSERT(i < m_num_species && j < m_num_species, "Out of bounds");
      return i + m_num_species * j;
    }

    std::size_t sym_index(std::size_t i, std::size_t j) const {
      return index(std::max(i, j), std::min(i, j));
    }
  };

}  // namespace otf
