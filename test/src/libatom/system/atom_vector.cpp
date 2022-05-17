#include "libatom/system/atom_vector.hpp"

#include <cmath>
#include <random>

#include "doctest/doctest.h"
#include "libatom/system/member.hpp"
#include "libatom/utils.hpp"

TEST_CASE("Atom vec") {
  using namespace otf;

  AtomVector<Position, AtomicNum> atoms;

  REQUIRE(atoms.size() == 0);

  ////// Add atom and test access

  atoms.emplace_back(Vec3<floating>::Zero(), 0);

  REQUIRE(atoms.size() == 1);

  REQUIRE(norm(atoms(Position{}, 0)) == 0);
  REQUIRE(atoms(AtomicNum{}, 0) == 0);

  atoms(Position{}) += 1;

  REQUIRE(std::abs(norm(atoms(Position{}, 0)) - std::sqrt(spatial_dims)) < 0.001);

  ////// Add second

  atoms.emplace_back(Vec3<floating>::Zero(), 0);

  REQUIRE(atoms.size() == 2);
  REQUIRE(atoms(AtomicNum{}, 1) == 0);

  ////// Add third

  atoms.emplace_back(Vec3<floating>::Zero(), 1);

  REQUIRE(atoms.size() == 3);
  REQUIRE(atoms(AtomicNum{}, 2) == 1);

  ////// Remove atoms

  atoms.shrink(2);

  REQUIRE(atoms.size() == 2);
  REQUIRE(atoms(AtomicNum{}).size() == 2);

  //// Re add third

  atoms.emplace_back(Vec3<floating>::Zero(), 5);

  REQUIRE(atoms.size() == 3);
  REQUIRE(atoms(AtomicNum{}, 2) == 5);
}