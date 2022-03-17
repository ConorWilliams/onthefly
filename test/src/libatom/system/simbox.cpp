#include "libatom/system/simbox.hpp"

#include <random>

#include "doctest/doctest.h"
#include "libatom/utils.hpp"

TEST_CASE("Simbox cannon_image") {
  using namespace otf;

  std::mt19937 gen(33);
  std::uniform_real_distribution<otf::float_t> dis(0, 1);

  auto vrand = [&] { return Vec<otf::float_t>::NullaryExpr([&]() { return dis(gen); }); };

  for (size_t i = 0; i < 100'000; i++) {
    //
    Vec<otf::float_t> extents = vrand() + 1;
    Vec<bool> periodic = vrand() < .5;

    OrthoSimBox box{extents, periodic};

    // Random points inside simbox
    Vec<otf::float_t> a = vrand() * extents;
    Vec<otf::float_t> b = vrand() * extents;

    // Displace by integral number of random extents in each periodic direction
    Vec<otf::float_t> b_prime = periodic.select(b + (10 * vrand()).floor() * extents, b);

    // Check in same position
    REQUIRE(norm(box.canon_image(b_prime) - a) - otf::norm(a - b) < 0.001);
  }
}

TEST_CASE("Simbox min_image") {}