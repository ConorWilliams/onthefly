#include "libatom/system/ortho_sim_box.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <random>

#include "libatom/utils.hpp"

TEST_CASE("Simbox cannon_image") {
  using namespace otf;

  std::mt19937 gen(33);
  std::uniform_real_distribution<floating> dis(0, 1);

  auto vrand = [&] { return Vec3<floating>::NullaryExpr([&]() { return dis(gen); }); };

  for (size_t i = 0; i < 100'000; i++) {
    //
    Vec3<floating> extents = vrand() + 1;
    Vec3<bool> periodic = vrand() < .5;

    OrthoSimBox box{extents, periodic};

    // Random points inside simbox
    Vec3<floating> a = vrand() * extents;
    Vec3<floating> b = vrand() * extents;

    // Displace by integral number of random extents in each periodic periodic direction
    Vec3<floating> b_prime = periodic.select(b + (10 * vrand()).floor() * extents, b);

    // Check in same position
    REQUIRE(std::abs(norm(box.canon_image(b_prime) - a) - otf::norm(a - b)) < 0.001);
  }
}

TEST_CASE("Simbox min_image") {
  //
  using namespace otf;

  OrthoSimBox box{Vec3<floating>::Constant(10), Vec3<bool>::Constant(true)};

  Vec3<floating> a = Vec3<floating>::Constant(1);
  Vec3<floating> b = Vec3<floating>::Constant(9);

  Vec3<floating> m = box.min_image(a, b);

  Vec3<floating> x = Vec3<floating>::Constant(-2);

  REQUIRE(norm(m - x) < 0.001);
}