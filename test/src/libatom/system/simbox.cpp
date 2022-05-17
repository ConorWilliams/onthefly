// #include "libatom/system/simbox.hpp"

// #include <cmath>
// #include <random>

// #include "doctest/doctest.h"
// #include "libatom/utils.hpp"

// TEST_CASE("Simbox cannon_image") {
//   using namespace otf;

//   std::mt19937 gen(33);
//   std::uniform_real_distribution<floating> dis(0, 1);

//   auto vrand = [&] { return Vec<floating>::NullaryExpr([&]() { return dis(gen); }); };

//   for (size_t i = 0; i < 100'000; i++) {
//     //
//     Vec<floating> extents = vrand() + 1;
//     Vec<bool> periodic = vrand() < .5;

//     OrthoSimBox box{extents, periodic};

//     // Random points inside simbox
//     Vec<floating> a = vrand() * extents;
//     Vec<floating> b = vrand() * extents;

//     // Displace by integral number of random extents in each periodic periodic direction
//     Vec<floating> b_prime = periodic.select(b + (10 * vrand()).floor() * extents, b);

//     // Check in same position
//     REQUIRE(std::abs(norm(box.canon_image(b_prime) - a) - otf::norm(a - b)) < 0.001);
//   }
// }

// TEST_CASE("Simbox min_image") {
//   //
//   using namespace otf;

//   OrthoSimBox box{Vec<floating>::Constant(10), Vec<bool>::Constant(true)};

//   Vec<floating> a = Vec<floating>::Constant(1);
//   Vec<floating> b = Vec<floating>::Constant(9);

//   Vec<floating> m = box.min_image(a, b);

//   Vec<floating> x = Vec<floating>::Constant(-2);

//   REQUIRE(norm(m - x) < 0.001);
// }