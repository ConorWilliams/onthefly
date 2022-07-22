#include <fmt/core.h>
#include <fmt/os.h>

#include <chrono>
#include <cmath>
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

SimCell init_cell(int n, floating l) {
  //   432
  //   It=-1 Lattice="17.1598406 0.0 0.0 0.0 17.1598406 0.0 0.0 0.0 17.1598406"

  SimCell cell({{l, l, l}, {true, true, true}});

  cell.destructive_resize(n);

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

  auto cell = init_cell(430, 17.1598406);

  std::ifstream file("/home/cdt1902/phd/P2021/data/v3h1/300k/olkmc.xyz", std::ios::in);

  env::EnvCell envs({5.2}, cell);

  for (std::size_t i = 0; i < 1; i++) {
    io::stream_xyz(file, cell);
  }

  envs.rebuild(cell, omp_get_max_threads());

  floating minx = 100;

  for (std::size_t i = 0; i < envs.size(); i++) {
    minx = std::min(minx, envs[i].fingerprint().r_min());
  }

  fmt::print("rmin={}\n", minx);

  floating start = 1e-14;
  floating stop = 1.0;

  std::size_t N = 1000;

  floating log_mult = std::log(stop / start) / (N - 1);

  auto fout = fmt::output_file("/home/cdt1902/phd/P2021/data/v3h1_envs.txt");

  fout.print("delta num_envs");

  for (std::size_t i = 0; i < N; i++) {
    floating delta = start * std::exp(log_mult * i);

    env::Catalogue cat({delta, 1.0});

    time_call("canon", [&] {
      for (std::size_t j = 0; j < cell.size(); j++) {
        // time_call("canon_find", &env::Catalogue::canon_find, cat, envs[j])

        if (!cat.canon_find(envs[j])) {
          cat.insert(envs[j]);
        }
      }
    });

    fout.print("\n{} {}", delta, cat.size());
    fout.flush();

    fmt::print("iter={} tot={} delta={}, keys={}\n", i, cat.size(), delta, cat.num_keys());
  }

  fmt::print("Done\n");

  return 0;
}