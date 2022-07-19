#include <fmt/core.h>
#include <fmt/os.h>

#include <cstddef>
#include <fstream>
#include <optional>

#include "libatom/atom.hpp"
#include "libatom/env/catalogue.hpp"
#include "libatom/env/topology.hpp"
#include "libatom/io/xyz.hpp"
#include "libatom/sim_cell.hpp"
#include "libatom/utils.hpp"

using namespace otf;

SimCell init_cell() {
  //   432
  //   It=-1 Lattice="17.1598406 0.0 0.0 0.0 17.1598406 0.0 0.0 0.0 17.1598406"

  SimCell cell({{17.1598406, 17.1598406, 17.1598406}, {true, true, true}});

  cell.destructive_resize(431);

  for (std::size_t i = 0; i < cell.size(); i++) {
    cell(Frozen{}, i) = false;
    cell(AtomicNum{}, i) = 26;
  }

  cell(AtomicNum{}, cell.size() - 1) = 1;  // Set last is hydrogen

  return cell;
}

using namespace otf;

auto main(int, char**) -> int {
  //

  fmt::print(stderr, "num_threads={}\n", omp_get_max_threads());

  auto cell = init_cell();

  std::ifstream file("/home/cdt1902/phd/P2021/data/v2h1/300k/olkmc.xyz", std::ios::in);

  int iter = 0;

  env::EnvCell envs({5.2}, cell);

  floating delta_max = 0.001;

  env::Catalogue cat({delta_max});

  auto fout = fmt::output_file("test.xyz");

  fmt::print("delta_max={}, num_atoms={}\n", delta_max, cell.size());

  while (!file.eof() && iter++ < 18) {
    //
    io::stream_xyz(file, cell);

    envs.rebuild(cell, omp_get_max_threads());

    std::size_t new_envs = 0;

    for (std::size_t i = 0; i < cell.size(); i++) {
      std::optional p = cat.canon_find(envs.env(i));

      if (!p) {
        cat.insert(envs.env(i));
        ++new_envs;
      }
    }

    fmt::print("{} {} {}\n", iter, new_envs, cat.size());
  }

  fmt::print("Done\n");

  return 0;
}