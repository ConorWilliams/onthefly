
#include "libatom/env/topology.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/env/geometry.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/utils.hpp"

namespace otf::env {

  void LocalEnv::rebuild(SimCell const &atoms, neighbour::List const &nl, std::size_t idx) {
    // Reuses our storage
    this->clear();
    m_key.clear();
    m_fingerprint.m_r_0j.clear();

    // Add central atom
    std::size_t centre_col = atoms.index_to_colour(idx);
    this->emplace_back({0, 0, 0}, centre_col, idx);
    m_key.emplace_back(centre_col, 1);

    nl.for_neighbours(idx, [&](std::size_t np, Vec3<double> const &dr) {
      //
      std::size_t n = nl.image_to_real(np);
      std::size_t n_col = atoms.index_to_colour(n);

      this->emplace_back(dr, n_col, n);

      m_fingerprint.m_r_0j.push_back(norm(dr));

      for (std::size_t i = 1; i < m_key.size(); i++) {
        if (m_key[i].first == n_col) {
          m_key[i].second++;
          return;
        }
      }
      // If here then we did not find an entry (no return) in m_key so we must create one.
      m_key.emplace_back(n_col, 1);
    });

    // Sort the colour orders
    std::sort(std::next(m_key.begin()), m_key.end());

    // Sort the r_0j part of the fingerprint
    std::sort(m_fingerprint.m_r_0j.begin(), m_fingerprint.m_r_0j.end());

    {  // Build r_ij part of the fingerprint

      m_fingerprint.m_r_ij.clear();

      for (std::size_t i = 1; i < size(); i++) {
        for (std::size_t j = 1; j < i; j++) {
          m_fingerprint.m_r_ij.push_back(norm((*this)[i](Position{}) - (*this)[j](Position{})));
        }
      }

      std::sort(m_fingerprint.m_r_ij.begin(), m_fingerprint.m_r_ij.end());
    }

    {  // Set COM == 0,0,0
      Vec3<floating> shift = com(*this) / (floating)size();

      for (auto &&elem : *this) {
        elem(Position{}) -= shift;
      }
    }
  }

  void EnvCell::rebuild(SimCell const &atoms, std::size_t num_threads) {
    // Always rebuild, no idea what has happened to atoms.
    m_nl.rebuild(atoms, num_threads);

    m_envs.resize(atoms.size());

#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (std::size_t i = 0; i < atoms.size(); i++) {
      m_envs[i].rebuild(atoms, m_nl, i);
    }
  }

}  // namespace otf::env