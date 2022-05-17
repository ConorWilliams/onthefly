#include "libatom/system/atomvector.hpp"

#include <cmath>
#include <random>

#include "doctest/doctest.h"
#include "libatom/utils.hpp"

TEST_CASE("Atom vec") {
  using namespace otf;

  AtomVector<Pos, AtomicNum> atoms;

  REQUIRE(atoms.size() == 0);

  //////

  atoms.emplace_back(Vec3<floating>::Zero(), 0);

  REQUIRE(atoms.size() == 1);

  REQUIRE(norm(atoms(Pos{}, 0)) == 0);
  REQUIRE(atoms(AtomicNum{}, 0) == 0);

  atoms(Pos{}) += 1;

  REQUIRE(std::abs(norm(atoms(Pos{}, 0)) - std::sqrt(spatial_dims)) < 0.001);

  //////

  atoms.emplace_back(Vec3<floating>::Zero(), 0);

  REQUIRE(atoms.size() == 2);
  REQUIRE(atoms(AtomicNum{}, 1) == 0);

  //////

  atoms.emplace_back(Vec3<floating>::Zero(), 1);

  REQUIRE(atoms.size() == 3);
  REQUIRE(atoms(AtomicNum{}, 2) == 1);
}