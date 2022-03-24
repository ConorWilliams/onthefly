#include <iostream>

#include "fmt/os.h"
#include "libatom/io/xyz.hpp"
#include "libatom/system/atomvector.hpp"
#include "libatom/system/simcell.hpp"

auto main(int, char **) -> int {
  //

  otf::SimCell cell({{10, 10, 10}, {true, true, false}});

  cell.active().emplace_back({1, 2, 3}, otf::Symbol{"Fe"});
  cell.active().emplace_back({5, 5, 5}, otf::Symbol{"H"});

  auto out = fmt::output_file("test.xyz", fmt::file::WRONLY | fmt::file::CREATE);

  otf::to_xyz(out, cell, 2.2);

  std::cout << "working\n";

  return 0;
}