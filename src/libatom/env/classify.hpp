#pragma once

#include <cstddef>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/atom.hpp"
#include "libatom/env/geometry.hpp"
#include "libatom/env/topology.hpp"
#include "libatom/neighbour/list.hpp"
#include "libatom/utils.hpp"

namespace otf::env {
  //

  class EnvCell {
  public:
    struct Options {
      // Radius of environment
      floating r_env;
    };

    EnvCell(Options const& opt, OrthoSimBox const& box) : m_nl(box, opt.r_env){};

    LocalEnv const& env(std::size_t i) const noexcept { return m_envs[i]; }

    LocalEnv& env(std::size_t i) noexcept { return m_envs[i]; }

    void rebuild(SimCell const& atoms, std::size_t num_threads);

  private:
    std::vector<LocalEnv> m_envs;

    neighbour::List m_nl;
  };

}  // namespace otf::env