#include "libatom/system/atom_array.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <random>

#include "libatom/system/member.hpp"
#include "libatom/utils.hpp"

TEST_CASE("Atom vec") {
  using namespace otf;

  AtomArray<Position, AtomicNum> atoms;

  REQUIRE(atoms.size() == 0);

  ////// Add atom and test access

  atoms.resize(1);

  atoms(Position{}, 0) = Vec3<floating>{0, 0, 0};
  atoms(AtomicNum{}, 0) = 0;

  REQUIRE(atoms.size() == 1);
  REQUIRE(norm(atoms(Position{}, 0)) == 0);
  REQUIRE(atoms(AtomicNum{}, 0) == 0);

  atoms(Position{}) += 1;

  REQUIRE(std::abs(norm(atoms(Position{}, 0)) - std::sqrt(spatial_dims)) < 0.001);

  ////// Add second

  atoms.resize(2);

  REQUIRE(atoms.size() == 2);
  REQUIRE(std::abs(norm(atoms(Position{}, 0)) - std::sqrt(spatial_dims)) < 0.001);
  REQUIRE(atoms(AtomicNum{}, 0) == 0);

  atoms(Position{}, 1) = Vec3<floating>{0, 0, 0};
  atoms(AtomicNum{}, 1) = 1;

  REQUIRE(atoms(AtomicNum{}, 1) == 1);
}