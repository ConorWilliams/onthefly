#include <fmt/core.h>
#include <fmt/os.h>

#include <array>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <optional>
#include <string>

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

  cell.destructive_resize(432);

  for (std::size_t i = 0; i < cell.size(); i++) {
    cell(Frozen{}, i) = false;
    cell(AtomicNum{}, i) = 26;
  }

  cell(AtomicNum{}, cell.size() - 1) = 1;  // Set last is hydrogen

  return cell;
}

using namespace otf;

auto main(int, char **) -> int {
  //

  fmt::print(stderr, "num_threads={}\n", omp_get_max_threads());

  auto cell = init_cell();

  std::ifstream file("/home/cdt1902/phd/P2021/data/v1h1/300k/olkmc.xyz", std::ios::in);

  floating D = 1e-1;

  env::EnvCell envs({5.2}, cell);

  env::Catalogue cat({D, 1.0});

  auto fout = fmt::output_file("/home/cdt1902/phd/P2021/data/fuzz_envs.1e-1.txt");

  fout.print("delta_max={}, num_atoms={}", D, cell.size());

  int iter = 0;

  while (!file.eof() && iter++ < 1'000) {
    //
    io::stream_xyz(file, cell);

    envs.rebuild(cell, omp_get_max_threads());

    for (std::size_t j = 0; j < cell.size(); j++) {
      if (!cat.canon_find(envs[j])) {
        cat.insert(envs[j]);
      }
    }

    fout.print("\n{} {}", iter, cat.size());

    fmt::print("iter={} tot={} keys={}\n", iter, cat.size(), cat.num_keys());

    fout.flush();

    if (iter % 100 == 0) {
      cat.optimize();
    }
  }

  fmt::print("Done\n");

  return 0;
}