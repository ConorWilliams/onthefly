#include "libatom/potentials/EAM/data.hpp"

#include <fstream>
#include <string_view>

#include "doctest/doctest.h"
#include "fmt/core.h"

using namespace otf;

TEST_CASE("Load DataEAM") {
  //

  if (std::ifstream file{"../data/wen.eam.fs"}; file.good()) {
    DataEAM data{file};
  } else if (std::ifstream file{"data/wen.eam.fs"}; file.good()) {
    DataEAM data{file};
  } else {
    REQUIRE(false);  // Could not locate data file.
  }
}
