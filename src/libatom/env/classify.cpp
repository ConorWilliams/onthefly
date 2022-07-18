

#include "libatom/env/classify.hpp"

#include <cstddef>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/env/geometry.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/utils.hpp"

namespace otf::env {

  void EnvCell::rebuild(SimCell const& atoms, std::size_t num_threads) {
    // Always rebuild, no idea what has happened to atoms.
    m_nl.rebuild(atoms, num_threads);

    m_envs.resize(atoms.size());

#pragma omp parallel for num_threads(num_threads) schedule(static)
    for (std::size_t i = 0; i < atoms.size(); i++) {
      m_envs[i].rebuild(atoms, m_nl, i);
    }
  }

}  // namespace otf::env