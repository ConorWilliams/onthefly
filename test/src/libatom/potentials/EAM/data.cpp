#include "libatom/potentials/EAM/data.hpp"

#include <doctest/doctest.h>
#include <fmt/core.h>

#include <fstream>

using namespace otf;

TEST_CASE("Load DataEAM") {
  //

  if (std::ifstream file{"../data/wen.eam.fs"}; file.good()) {
    DataEAM data{std::move(file)};
  } else if (std::ifstream file{"data/wen.eam.fs"}; file.good()) {
    DataEAM data{std::move(file)};
  } else {
    fmt::print(stderr, "Could not locate EAM file\n");
  }
}
