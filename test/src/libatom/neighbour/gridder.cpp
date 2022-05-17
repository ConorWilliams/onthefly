// #include "libatom/neighbour/gridder.hpp"

// #include <array>
// #include <cmath>
// #include <random>

// #include "doctest/doctest.h"
// #include "libatom/system/simbox.hpp"
// #include "libatom/utils.hpp"

// TEST_CASE("Gridder") {
//   otf::OrthoSimBox box{{10, 10, 10}, {true, true, true}};

//   otf::Gridder grid;

//   grid.compute_neigh_cells(box, 3);

//   REQUIRE(grid.cell_idx(otf::Vec<otf::floating>{0, 0, 0} + grid.origin()) == 1 + 1 * 5 + 1 * 5 *
//   5); REQUIRE(grid.cell_idx(otf::Vec<otf::floating>{0, 0, 0} + grid.origin()) == 1 + 1 * 5 + 1 *
//   5 * 5);

//   REQUIRE(grid.cell_idx(otf::Vec<otf::floating>{5, 5, 5} + grid.origin()) == 2 + 2 * 5 + 2 * 5 *
//   5);

//   REQUIRE(grid.cell_idx(otf::Vec<otf::floating>{9.999, 9.999, 9.999} + grid.origin())
//           == 3 + 3 * 5 + 3 * 5 * 5);

//   std::array neigh = grid.neigh_cells(1 + 1 * 5 + 1 * 5 * 5);

//   REQUIRE(neigh[0] == (1 - 1) + (1 - 1) * 5 + (1 - 1) * 5 * 5);
//   REQUIRE(neigh[1] == (1 - 0) + (1 - 1) * 5 + (1 - 1) * 5 * 5);
//   REQUIRE(neigh[2] == (1 + 1) + (1 - 1) * 5 + (1 - 1) * 5 * 5);

//   REQUIRE(neigh[3] == (1 - 1) + (1 - 0) * 5 + (1 - 1) * 5 * 5);
//   REQUIRE(neigh[4] == (1 - 0) + (1 - 0) * 5 + (1 - 1) * 5 * 5);
//   REQUIRE(neigh[5] == (1 + 1) + (1 - 0) * 5 + (1 - 1) * 5 * 5);

//   REQUIRE(neigh[6] == (1 - 1) + (1 + 1) * 5 + (1 - 1) * 5 * 5);
//   REQUIRE(neigh[7] == (1 - 0) + (1 + 1) * 5 + (1 - 1) * 5 * 5);
//   REQUIRE(neigh[8] == (1 + 1) + (1 + 1) * 5 + (1 - 1) * 5 * 5);

//   REQUIRE(neigh[9] == (1 - 1) + (1 - 1) * 5 + (1 - 0) * 5 * 5);
//   REQUIRE(neigh[10] == (1 - 0) + (1 - 1) * 5 + (1 - 0) * 5 * 5);
//   REQUIRE(neigh[11] == (1 + 1) + (1 - 1) * 5 + (1 - 0) * 5 * 5);

//   REQUIRE(neigh[12] == (1 - 1) + (1 - 0) * 5 + (1 - 0) * 5 * 5);
//   //   REQUIRE(neigh[13] == (1 - 0) + (1 - 0) * 5 + (1 - 0) * 5 * 5);
//   REQUIRE(neigh[13] == (1 + 1) + (1 - 0) * 5 + (1 - 0) * 5 * 5);

//   REQUIRE(neigh[14] == (1 - 1) + (1 + 1) * 5 + (1 - 0) * 5 * 5);
//   REQUIRE(neigh[15] == (1 - 0) + (1 + 1) * 5 + (1 - 0) * 5 * 5);
//   REQUIRE(neigh[16] == (1 + 1) + (1 + 1) * 5 + (1 - 0) * 5 * 5);

//   REQUIRE(neigh[17] == (1 - 1) + (1 - 1) * 5 + (1 + 1) * 5 * 5);
//   REQUIRE(neigh[18] == (1 - 0) + (1 - 1) * 5 + (1 + 1) * 5 * 5);
//   REQUIRE(neigh[19] == (1 + 1) + (1 - 1) * 5 + (1 + 1) * 5 * 5);

//   REQUIRE(neigh[20] == (1 - 1) + (1 - 0) * 5 + (1 + 1) * 5 * 5);
//   REQUIRE(neigh[21] == (1 - 0) + (1 - 0) * 5 + (1 + 1) * 5 * 5);
//   REQUIRE(neigh[22] == (1 + 1) + (1 - 0) * 5 + (1 + 1) * 5 * 5);

//   REQUIRE(neigh[23] == (1 - 1) + (1 + 1) * 5 + (1 + 1) * 5 * 5);
//   REQUIRE(neigh[24] == (1 - 0) + (1 + 1) * 5 + (1 + 1) * 5 * 5);
//   REQUIRE(neigh[25] == (1 + 1) + (1 + 1) * 5 + (1 + 1) * 5 * 5);
// }