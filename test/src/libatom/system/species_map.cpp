#include "libatom/system/species_map.hpp"

#include <cmath>
#include <random>

#include "doctest/doctest.h"
#include "libatom/utils.hpp"

TEST_CASE("SpeciesMap") {
  using namespace otf;

  SpeciesMap map;

  REQUIRE(map.size() == 0);

  //////

  REQUIRE(!map.species2z("not in"));

  auto a = map.species2z_or_insert("a");
  auto b = map.species2z_or_insert("b");

  REQUIRE(a == 0);
  REQUIRE(b == 1);

  REQUIRE(map.species2z("a"));
  REQUIRE(map.species2z("b"));

  REQUIRE(*map.species2z("a") == a);
  REQUIRE(*map.species2z("b") == b);

  REQUIRE(!map.species2z("not in"));

  REQUIRE(map.z2species(a) == "a");
  REQUIRE(map.z2species(b) == "b");
}