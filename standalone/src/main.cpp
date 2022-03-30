

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

  vec.active().emplace_back({0, 0, 1}, otf::Symbol{"H"});

  {
    std::fstream s{"dump.bin", s.binary | s.trunc | s.out};

    // we cannot use quick serialization function, because streams cannot use writtenBytesCount
    //   method

    otf::dump_binary(s, vec);

    vec.active().x().col(0) += 1;

    otf::dump_binary(s, vec);
  }

  STACK();

  {
    std::fstream s("dump.bin", s.binary | s.in);

    otf::SimCell res{{{1, 9, 1}, {false, false, false}}};

    while (true) {
      bool last = stream_binary(s, res);

      // Do something with res

      fmt::print("{}\n", fmt::join(res.active().x().col(0), " "));

      if (last) {
        break;
      }
    }

    // }
  }

  //   int a = 1;

  otf::timeit("tdest", []() {});

  //

  std::cout << "working\n";

  return 0;
}