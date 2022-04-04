#include "libatom/system/atomvector.hpp"

#include <cmath>
#include <random>

#include "doctest/doctest.h"
#include "libatom/utils.hpp"

TEST_CASE("Atom vec") {
  using namespace otf;

  AtomVector atoms;

  REQUIRE(atoms.size() == 0);

  //////

  atoms.emplace_back(Vec<flt_t>::Zero(), 0);

  REQUIRE(atoms.size() == 1);

  REQUIRE(norm(atoms.x().col(0)) == 0);
  REQUIRE(atoms.z()[0] == 0);

  atoms.x() += 1;

  REQUIRE(std::abs(norm(atoms.x().col(0)) - std::sqrt(spatial_dims)) < 0.001);

  //////

  atoms.emplace_back(Vec<flt_t>::Zero(), 0);

  REQUIRE(atoms.size() == 2);
  REQUIRE(atoms.z()[1] == 0);

  //////

  atoms.emplace_back(Vec<flt_t>::Zero(), 1);

  REQUIRE(atoms.size() == 3);
  REQUIRE(atoms.z()[2] == 1);
}