#include "libatom/system/simbox.hpp"

#include <cmath>
#include <random>

#include "doctest/doctest.h"
#include "libatom/utils.hpp"

TEST_CASE("Simbox cannon_image") {
  using namespace otf;

  std::mt19937 gen(33);
  std::uniform_real_distribution<flt_t> dis(0, 1);

  auto vrand = [&] { return Vec<flt_t>::NullaryExpr([&]() { return dis(gen); }); };

  for (size_t i = 0; i < 100'000; i++) {
    //
    Vec<flt_t> extents = vrand() + 1;
    Vec<bool> periodic = vrand() < .5;

    OrthoSimBox box{extents, periodic};

    // Random points inside simbox
    Vec<flt_t> a = vrand() * extents;
    Vec<flt_t> b = vrand() * extents;

    // Displace by integral number of random extents in each periodic periodic direction
    Vec<flt_t> b_prime = periodic.select(b + (10 * vrand()).floor() * extents, b);

    // Check in same position
    REQUIRE(std::abs(norm(box.canon_image(b_prime) - a) - otf::norm(a - b)) < 0.001);
  }
}

TEST_CASE("Simbox min_image") {
  //
  using namespace otf;

  OrthoSimBox box{Vec<flt_t>::Constant(10), Vec<bool>::Constant(true)};

  Vec<flt_t> a = Vec<flt_t>::Constant(1);
  Vec<flt_t> b = Vec<flt_t>::Constant(9);

  Vec<flt_t> m = box.min_image(a, b);

  Vec<flt_t> x = Vec<flt_t>::Constant(-2);

  REQUIRE(norm(m - x) < 0.001);
}