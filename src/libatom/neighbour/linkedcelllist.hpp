
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "libatom/asserts.hpp"
#include "libatom/neighbour/gridder.hpp"
#include "libatom/system/simbox.hpp"
#include "libatom/system/simcell.hpp"
#include "libatom/utils.hpp"

namespace otf {

  /**
   * @brief  Abstraction to apply a function f as f(n, r_{an}, dr3) for every neighbour "n" of an
   atom
   * "a".
   *
   * Internally builds a vector of atoms, each with a Data object embedded within them, which
   contains
   * the active atoms, boundary atoms and any ghost atoms required if the supercell is periodic.
   *
   * @tparam Data
   */

  template <typename T> class LinkedCellList {
  public:
    /**
     * @brief The atom type used by LinkedCellList, can be customised through template parameter
     *
     * After constructing a LinkedCellList you will be able to iterate over all pairs of neighbours
     * within some cut-off, the neighbors will be periodically resolved (some active, some frozen
     * and some "ghosts") and will be of this type.
     *
     * @tparam T Atom will derive from T publicly and can be used to inject members into Atom
     */
    class Atom : public T {
    public:
      /**
       * @brief The coordinates of the atom
       */
      Vec<float_t> x;

    private:
      static_assert(!std::is_fundamental_v<T>);
      static_assert(std::is_trivially_destructible_v<T>);
      static_assert(std::is_default_constructible_v<T>);
      static_assert(std::is_assignable_v<T &, T &>);

      friend class LinkedCellList;

      Atom(Vec<double> const &x, T const &base) : T{base}, x{x}, next{nullptr} {}

      Atom *next;

      T &data() { return static_cast<T &>(*this); }
    };

    /**
     * @brief Build the linked cell list, binning atoms into orthogananl cells (min "rcut" in each
     * dimension).
     *
     * @param Fn A function that must have the signature fn(SimCell const &, std::size_t i) -> T,
     * used to construct the injected T member of each atom.
     */
    template <typename F> void construct(SimCell const &cell, double rcut, F &&fn) {
      STACK();
      m_active = cell.active.size();
      m_size = cell.size();

      grid.compute_neigh_cells(cell.box, rcut);

      static_assert(std::is_trivially_destructible_v<Atom>, "");
      m_list.clear();  // Should be order one by ^

      // Copy active atoms into list IN ORDER
      for (std::size_t idx = 0; idx < cell.active.size(); ++idx) {
        m_list.emplace_back(cell.box.canon_image(cell.active.x().col(idx)) + grid.origin(),
                            fn(cell, idx));
      }

      // Copy frozen atoms into list IN ORDER
      for (std::size_t idx = 0; idx < cell.frozen.size(); ++idx) {
        m_list.emplace_back(cell.box.canon_image(cell.frozen.x().col(idx)) + grid.origin(),
                            fn(cell, idx));
      }

      m_list_idx.clear();  // Should be order one

      for (auto &&atom : m_list) {
        auto idx = grid.

                   m_list_idx.emplace_back(grid.cell_idx(const Vec<flt_t> &x))
      }

      make_ghosts();

      // Update head;
      _head.assign(_prod_shape[2], nullptr);

      for (auto &atom : _list) {
        atom.next = std::exchange(_head[lambda(atom)], &atom);
      }
    }

    void broadcast_ghost_data();

    // Call f(i, r, dr) for each neighbour (index i) of idx
    // template <typename F> void neigh_reduce(neigh_atom *atom, F &&f);

  private:
    std::size_t m_active;
    std::size_t m_size;

    Gridder grid;

    std::vector<Atom *> m_head;
    std::vector<Atom> m_list;
    std::vector<int> m_list_idx;

    // // Rapaport p.18
    // void make_ghosts();
  };

  //   template <typename Data> void LinkedCellList<Data>::load(double rcut, Supercell const &cell)
  //   {
  //     _active = cell.activ.size();
  //     _size = cell.size();

  //     load_simbox(rcut, cell);

  //     static_assert(std::is_trivially_destructible_v<neigh_atom>, "");
  //     _list.clear();  // Should be order one by ^

  //     // Copy active atoms into list IN ORDER
  //     for (std::size_t idx = 0; idx < cell.activ.size(); ++idx) {
  //       _list.emplace_back(cell.canonicle_image(cell.activ[idx].vec) + _cell,
  //       cell.activ[idx].col);
  //     }

  //     // Copy active atoms into list IN ORDER
  //     for (std::size_t idx = 0; idx < cell.bound.size(); ++idx) {
  //       _list.emplace_back(cell.canonicle_image(cell.bound[idx].vec) + _cell,
  //       cell.bound[idx].col);
  //     }

  // #if !defined(NDEBUG) && !defined(OLKMC_NO_DEBUG)
  //     // Verify all atoms mapped into Supercell correctly
  //     for (auto &&atom : _list) {
  //       int i = atom.vec[0] * _inv_cell[0];
  //       int j = atom.vec[1] * _inv_cell[1];
  //       int k = atom.vec[2] * _inv_cell[2];

  //       CHECK(i >= 1, "x-bound in canonicalize");
  //       CHECK(j >= 1, "y-bound in canonicalize");
  //       CHECK(k >= 1, "z-bound in canonicalize");

  //       CHECK(i < _shape[0] - 1, "x-bound in canonicalize");
  //       CHECK(j < _shape[1] - 1, "y-bound in canonicalize");
  //       CHECK(k < _shape[2] - 1, "z-bound in canonicalize");
  //     }
  // #endif

  //     make_ghosts();

  //     // Update head;
  //     _head.assign(_prod_shape[2], nullptr);

  //     for (auto &atom : _list) {
  //       atom.next = std::exchange(_head[lambda(atom)], &atom);
  //     }
  //   }

  //   template <typename Data> void LinkedCellList<Data>::broadcast_ghost_data() {
  //     std::size_t ghost_count = 0;

  //     for (int i = 0; i < 3; ++i) {
  //       if (_box.periodic & (1U << i)) {
  //         int const end = _size + ghost_count;

  //         for (int j = 0; j < end; ++j) {
  //           if (_list[j].vec[i] < _cell[i] + _rcut) {
  //             _list[_size + ghost_count++].data() = _list[j].data();
  //           }

  //           if (_list[j].vec[i] >= _box.extents[i] + _cell[i] - _rcut) {
  //             _list[_size + ghost_count++].data() = _list[j].data();
  //           }
  //         }
  //       }
  //     }
  //   }

  //   // Call f(i, r, dr) for each neighbour (index i) of idx
  //   template <typename Data> template <typename F>
  //   void LinkedCellList<Data>::neigh_reduce(neigh_atom *atom, F &&f) {
  //     long const lam = lambda(*atom);

  //     neigh_atom *neigh = _head[lam];

  //     // In same cell must check not-self
  //     while (neigh != nullptr) {
  //       if (neigh != atom) {
  //         Vec3<double> dr = atom->vec - neigh->vec;

  //         double r_sq = norm_sq(dr);

  //         if (r_sq < _rcut_sq) {
  //           if constexpr (std::is_invocable_v<F, neigh_atom *, double, Vec3<double>>) {
  //             f(neigh, std::sqrt(r_sq), dr);
  //           } else if constexpr (std::is_invocable_v<F, neigh_atom *, double>) {
  //             f(neigh, std::sqrt(r_sq));
  //           } else {
  //             f(neigh);
  //           }
  //         }
  //       }
  //       neigh = neigh->next;
  //     }

  //     // In adjacent cells -- don't need to check against self
  //     for (auto off : _neigh_stride) {
  //       neigh = _head[lam + off];

  //       while (neigh != nullptr) {
  //         Vec3<double> dr = atom->vec - neigh->vec;

  //         double r_sq = norm_sq(dr);

  //         if (r_sq < _rcut_sq) {
  //           if constexpr (std::is_invocable_v<F, neigh_atom *, double, Vec3<double>>) {
  //             f(neigh, std::sqrt(r_sq), dr);
  //           } else if constexpr (std::is_invocable_v<F, neigh_atom *, double>) {
  //             f(neigh, std::sqrt(r_sq));
  //           } else {
  //             f(neigh);
  //           }
  //         }

  //         neigh = neigh->next;
  //       }
  //     }

  //     return;
  //   }

  //   template <typename Data> void LinkedCellList<Data>::load_simbox(double rcut, Simbox const
  //   &box) {
  //     if (_rcut == rcut && box == _box) {
  //       return;
  //     }

  //     _rcut = rcut;
  //     _rcut_sq = _rcut * _rcut;

  //     _box = box;

  //     _shape = 2 + (_box.extents / _rcut).cast<int>();
  //     _cell = _box.extents / (_box.extents / rcut).floor();
  //     _inv_cell = 1.0 / _cell;

  //     // Sanity checks
  //     CHECK(_rcut > 0, "rcut is negative");
  //     CHECK((box.extents >= _rcut).all(), "rcut is too big");

  //     // Cumprod _shape
  //     _prod_shape[0] = _shape[0];
  //     _prod_shape[1] = _shape[0] * _shape[1];
  //     _prod_shape[2] = _shape[0] * _shape[1] * _shape[2];

  //     // Compute neighbour stride offsets
  //     int idx = 0;
  //     for (auto k : {-1, 0, 1}) {
  //       for (auto j : {-1, 0, 1}) {
  //         for (auto i : {-1, 0, 1}) {
  //           if (i != 0 || j != 0 || k != 0) {
  //             _neigh_stride[idx++] = i + j * _shape[0] + k * _prod_shape[1];
  //           }
  //         }
  //       }
  //     }
  //   }

  //   // Makes ghost atoms, Rapaport p.18
  //   template <typename Data> void LinkedCellList<Data>::make_ghosts() {
  //     for (int i = 0; i < 3; ++i) {
  //       // Only make ghosts if axis is periodic
  //       if (_box.periodic & (1U << i)) {
  //         int const end = _list.size();

  //         for (int j = 0; j < end; ++j) {
  //           // Must cache atom for push_back that resizes vector
  //           neigh_atom atom = _list[j];

  //           if (atom.vec[i] < _cell[i] + _rcut) {
  //             _list.push_back(atom);
  //             _list.back().vec[i] += _box.extents[i];

  //             CHECK(static_cast<int>(_list.back().vec[i] * _inv_cell[i]) == _shape[i] - 1,
  //                   "ghost maths error");
  //           }

  //           if (atom.vec[i] >= _box.extents[i] + _cell[i] - _rcut) {
  //             _list.push_back(atom);
  //             _list.back().vec[i] -= _box.extents[i];

  //             CHECK(static_cast<int>(_list.back().vec[i] * _inv_cell[i]) == 0, "ghost maths
  //             error");
  //           }
  //         }
  //       }
  //     }
  //   }

}  // namespace otf