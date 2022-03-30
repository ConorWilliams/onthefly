
// #pragma once

// #include <array>
// #include <cstddef>
// #include <type_traits>
// #include <utility>
// #include <variant>
// #include <vector>

// #include "config.hpp"
// #include "supercell.hpp"

// /**
//  * @brief  Abstraction to apply a function f as f(n, r_{an}, dr3) for every neighbour "n" of an
//  atom
//  * "a".
//  *
//  * Internally builds a vector of atoms, each with a Data object embedded within them, which
//  contains
//  * the active atoms, boundary atoms and any ghost atoms required if the supercell is periodic.
//  *
//  * @tparam Data
//  */
// template <typename Data> class NeighReduce {
// public:
//   struct neigh_atom : Data {
//     static_assert(!std::is_fundamental_v<Data>, "");
//     static_assert(std::is_trivially_destructible_v<Data>, "");
//     static_assert(std::is_default_constructible_v<Data>, "");
//     static_assert(std::is_assignable_v<Data &, Data &>, "");

//     Vec3<double> vec;
//     Colour col;
//     neigh_atom *next;

//     neigh_atom(Vec3<double> const &vec, Colour const &col) : vec{vec}, col{col}, next{nullptr} {}

//     Data &data() { return static_cast<Data &>(*this); }
//   };

//   //   neigh_atom *operator[](std::size_t i) { return _list.data() + i; }

//   neigh_atom *begin_activ() { return _list.data(); }
//   neigh_atom *begin_bound() { return _list.data() + _active; }
//   neigh_atom *begin_ghost() { return _list.data() + _size; }

//   void load(double rcut, Supercell const &cell);

//   void broadcast_ghost_data();

//   // Call f(i, r, dr) for each neighbour (index i) of idx
//   template <typename F> void neigh_reduce(neigh_atom *atom, F &&f);

// private:
//   std::size_t _active;
//   std::size_t _size;

//   Simbox _box;

//   double _rcut;
//   double _rcut_sq;

//   Vec3<int> _shape;
//   Vec3<int> _prod_shape;

//   Vec3<double> _cell;
//   Vec3<double> _inv_cell;

//   std::vector<neigh_atom *> _head;
//   std::vector<neigh_atom> _list;

//   std::array<int, 3 * 3 * 3 - 1> _neigh_stride;

//   // Maps canonical image position to integer
//   int lambda(neigh_atom const &atom) const;

//   void load_simbox(double rcut, Simbox const &box);

//   // Rapaport p.18
//   void make_ghosts();
// };
// template <typename Data> void NeighReduce<Data>::load(double rcut, Supercell const &cell) {
//   _active = cell.activ.size();
//   _size = cell.size();

//   load_simbox(rcut, cell);

//   static_assert(std::is_trivially_destructible_v<neigh_atom>, "");
//   _list.clear();  // Should be order one by ^

//   // Copy active atoms into list IN ORDER
//   for (std::size_t idx = 0; idx < cell.activ.size(); ++idx) {
//     _list.emplace_back(cell.canonicle_image(cell.activ[idx].vec) + _cell, cell.activ[idx].col);
//   }

//   // Copy active atoms into list IN ORDER
//   for (std::size_t idx = 0; idx < cell.bound.size(); ++idx) {
//     _list.emplace_back(cell.canonicle_image(cell.bound[idx].vec) + _cell, cell.bound[idx].col);
//   }

// #if !defined(NDEBUG) && !defined(OLKMC_NO_DEBUG)
//   // Verify all atoms mapped into Supercell correctly
//   for (auto &&atom : _list) {
//     int i = atom.vec[0] * _inv_cell[0];
//     int j = atom.vec[1] * _inv_cell[1];
//     int k = atom.vec[2] * _inv_cell[2];

//     CHECK(i >= 1, "x-bound in canonicalize");
//     CHECK(j >= 1, "y-bound in canonicalize");
//     CHECK(k >= 1, "z-bound in canonicalize");

//     CHECK(i < _shape[0] - 1, "x-bound in canonicalize");
//     CHECK(j < _shape[1] - 1, "y-bound in canonicalize");
//     CHECK(k < _shape[2] - 1, "z-bound in canonicalize");
//   }
// #endif

//   make_ghosts();

//   // Update head;
//   _head.assign(_prod_shape[2], nullptr);

//   for (auto &atom : _list) {
//     atom.next = std::exchange(_head[lambda(atom)], &atom);
//   }
// }

// template <typename Data> void NeighReduce<Data>::broadcast_ghost_data() {
//   std::size_t ghost_count = 0;

//   for (int i = 0; i < 3; ++i) {
//     if (_box.periodic & (1U << i)) {
//       int const end = _size + ghost_count;

//       for (int j = 0; j < end; ++j) {
//         if (_list[j].vec[i] < _cell[i] + _rcut) {
//           _list[_size + ghost_count++].data() = _list[j].data();
//         }

//         if (_list[j].vec[i] >= _box.extents[i] + _cell[i] - _rcut) {
//           _list[_size + ghost_count++].data() = _list[j].data();
//         }
//       }
//     }
//   }
// }

// // Call f(i, r, dr) for each neighbour (index i) of idx
// template <typename Data> template <typename F>
// void NeighReduce<Data>::neigh_reduce(neigh_atom *atom, F &&f) {
//   long const lam = lambda(*atom);

//   neigh_atom *neigh = _head[lam];

//   // In same cell must check not-self
//   while (neigh != nullptr) {
//     if (neigh != atom) {
//       Vec3<double> dr = atom->vec - neigh->vec;

//       double r_sq = norm_sq(dr);

//       if (r_sq < _rcut_sq) {
//         if constexpr (std::is_invocable_v<F, neigh_atom *, double, Vec3<double>>) {
//           f(neigh, std::sqrt(r_sq), dr);
//         } else if constexpr (std::is_invocable_v<F, neigh_atom *, double>) {
//           f(neigh, std::sqrt(r_sq));
//         } else {
//           f(neigh);
//         }
//       }
//     }
//     neigh = neigh->next;
//   }

//   // In adjacent cells -- don't need to check against self
//   for (auto off : _neigh_stride) {
//     neigh = _head[lam + off];

//     while (neigh != nullptr) {
//       Vec3<double> dr = atom->vec - neigh->vec;

//       double r_sq = norm_sq(dr);

//       if (r_sq < _rcut_sq) {
//         if constexpr (std::is_invocable_v<F, neigh_atom *, double, Vec3<double>>) {
//           f(neigh, std::sqrt(r_sq), dr);
//         } else if constexpr (std::is_invocable_v<F, neigh_atom *, double>) {
//           f(neigh, std::sqrt(r_sq));
//         } else {
//           f(neigh);
//         }
//       }

//       neigh = neigh->next;
//     }
//   }

//   return;
// }

// // Maps canonical image position to integer
// template <typename Data> int NeighReduce<Data>::lambda(neigh_atom const &atom) const {
//   int i = atom.vec[0] * _inv_cell[0];
//   int j = atom.vec[1] * _inv_cell[1];
//   int k = atom.vec[2] * _inv_cell[2];

//   CHECK(i >= 0, "x-bound");
//   CHECK(j >= 0, "y-bound");
//   CHECK(k >= 0, "z-bound");

//   CHECK(i < _shape[0], "x-bound");
//   CHECK(j < _shape[1], "y-bound");
//   CHECK(k < _shape[2], "z-bound");

//   return i + j * _prod_shape[0] + k * _prod_shape[1];
// }

// template <typename Data> void NeighReduce<Data>::load_simbox(double rcut, Simbox const &box) {
//   if (_rcut == rcut && box == _box) {
//     return;
//   }

//   _rcut = rcut;
//   _rcut_sq = _rcut * _rcut;

//   _box = box;

//   _shape = 2 + (_box.extents / _rcut).cast<int>();
//   _cell = _box.extents / (_box.extents / rcut).floor();
//   _inv_cell = 1.0 / _cell;

//   // Sanity checks
//   CHECK(_rcut > 0, "rcut is negative");
//   CHECK((box.extents >= _rcut).all(), "rcut is too big");

//   // Cumprod _shape
//   _prod_shape[0] = _shape[0];
//   _prod_shape[1] = _shape[0] * _shape[1];
//   _prod_shape[2] = _shape[0] * _shape[1] * _shape[2];

//   // Compute neighbour stride offsets
//   int idx = 0;
//   for (auto k : {-1, 0, 1}) {
//     for (auto j : {-1, 0, 1}) {
//       for (auto i : {-1, 0, 1}) {
//         if (i != 0 || j != 0 || k != 0) {
//           _neigh_stride[idx++] = i + j * _shape[0] + k * _prod_shape[1];
//         }
//       }
//     }
//   }
// }

// // Makes ghost atoms, Rapaport p.18
// template <typename Data> void NeighReduce<Data>::make_ghosts() {
//   for (int i = 0; i < 3; ++i) {
//     // Only make ghosts if axis is periodic
//     if (_box.periodic & (1U << i)) {
//       int const end = _list.size();

//       for (int j = 0; j < end; ++j) {
//         // Must cache atom for push_back that resizes vector
//         neigh_atom atom = _list[j];

//         if (atom.vec[i] < _cell[i] + _rcut) {
//           _list.push_back(atom);
//           _list.back().vec[i] += _box.extents[i];

//           CHECK(static_cast<int>(_list.back().vec[i] * _inv_cell[i]) == _shape[i] - 1,
//                 "ghost maths error");
//         }

//         if (atom.vec[i] >= _box.extents[i] + _cell[i] - _rcut) {
//           _list.push_back(atom);
//           _list.back().vec[i] -= _box.extents[i];

//           CHECK(static_cast<int>(_list.back().vec[i] * _inv_cell[i]) == 0, "ghost maths error");
//         }
//       }
//     }
//   }
// }