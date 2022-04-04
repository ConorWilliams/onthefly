

#include <fmt/core.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>

#include "fmt/ranges.h"
#include "libatom/io/binary.hpp"
#include "libatom/io/xyz.hpp"
#include "libatom/system/atomvector.hpp"
#include "libatom/system/simbox.hpp"
#include "libatom/system/simcell.hpp"
#include "libatom/utils.hpp"

auto main(int, char **) -> int {
  //

  otf::SimCell vec{{{1, 1, 1}, {false, false, false}}};

  auto fe_id = vec.map.species2z_or_insert("Fe");

  vec.active.emplace_back({0, 0, 1}, fe_id);

  {
    std::fstream s{"dump.bin", s.binary | s.trunc | s.out};

    otf::dump_binary(s, vec);

    vec.active.x().col(0) += 1;

    otf::dump_binary(s, vec);
  }

  fmt::print("{}\n", vec.map.z2species(fe_id));

  STACK();

  {
    std::fstream s("dump.bin", s.binary | s.in);

    otf::SimCell res{{{1, 9, 1}, {false, false, false}}};

    while (true) {
      bool last = stream_binary(s, res);

      // Do something with res

      fmt::print("{}\n", fmt::join(res.active.x().col(0), " "));

      VERIFY(res.map.species2z("Fe") && *res.map.species2z("Fe") == 0, "WOW");

      if (last) {
        break;
      }
    }

    // }
  }

  //   //   int a = 1;

  otf::timeit("tdest", []() {});

  //

  std::cout << "working\n";

  return 0;
}