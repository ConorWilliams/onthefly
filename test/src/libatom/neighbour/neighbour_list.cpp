#include "libatom/neighbour/neighbour_list.hpp"

#include "doctest/doctest.h"
#include "libatom/system/member.hpp"
#include "libatom/system/sim_cell.hpp"
#include "libatom/utils.hpp"

TEST_CASE("Neighbour list") {
  using namespace otf;

  // Set up sim cell

  SimCell atoms({{1, 1, 1}, {true, true, true}});

  atoms.resize(1);

  atoms(Position{}, 0) = Vec3<floating>{0.5, 0.5, 0.5};
  atoms(AtomicNum{}, 0) = 1;
  atoms(Frozen{}, 0) = false;

  NeighbourCell neigh;

  neigh.rebuild_neighbour_lists(atoms, 1);
}