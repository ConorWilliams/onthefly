#pragma once

#include <cstdint>
#include <iomanip>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "aligned_allocator.hpp"
#include "cereal/types/base_class.hpp"
#include "cereal/types/vector.hpp"
#include "config.hpp"
#include "toml++/toml.h"
#include "utility.hpp"

// Small struct to pack atomic_number and bound || active state into 32 bits
struct Colour {
    static constexpr std::size_t num_states = 2;

    enum [[deprecated]] : std::uint16_t{Fe = 0, H = 1};

    enum : std::uint16_t { activ = 0, bound = 1, vacant = 2 };

    std::uint16_t atomic;
    std::uint16_t state;

    // Implicitly combine atomic + state into colour
    constexpr operator std::size_t() const { return atomic * num_states + state; }

    // Returns maximum possible value that Colour could implicitly convert to
    static constexpr std::size_t max() { return num_states * NUM_ATOM_SPECIES; }

    template <class Archive> void serialize(Archive &ar) { ar(atomic, state); }
};

class AtomVector {
  private:
    struct AtomView {
        Eigen::Map<Vec3<double>> vec;
        Colour &col;
    };

    struct ConstAtomView {
        Eigen::Map<Vec3<double> const> vec;
        Colour const &col;
    };

  public:
    // Minimal std::vector API
    std::size_t size() const { return _col.size(); }

    // Fetch a proxy object that behaves like an atom
    ConstAtomView operator[](std::size_t n) const {
        CHECK(n < size(), "Bad access");
        return {Eigen::Map<Vec3<double> const>{_vec.data() + 3 * n}, _col[n]};
    }
    AtomView operator[](std::size_t n) {
        CHECK(n < size(), "Bad access");
        return {Eigen::Map<Vec3<double>>{_vec.data() + 3 * n}, _col[n]};
    }

    void emplace_back(Vec3<double> const &vec, Colour const &col) {
        _vec.insert(_vec.end(), vec.begin(), vec.end());
        _col.insert(_col.end(), col);
    }

    void clear() {
        _vec.clear();
        _col.clear();
    }

    // Eigen API

    // Fetch an Eigen view into the positions of the active atoms
    Eigen::Map<VecN<double>, OLKMC_EIGEN_ALIGN> view() {
        return {_vec.data(), static_cast<Eigen::Index>(_vec.size())};
    }

    Eigen::Map<VecN<double> const, OLKMC_EIGEN_ALIGN> view() const {
        return {_vec.data(), static_cast<Eigen::Index>(_vec.size())};
    }

    // For initialising VecN<double> large enough to represent all atoms
    VecN<double> zero_like() const { return VecN<double>::Zero(_vec.size()); }

    template <class Archive> void serialize(Archive &ar) { ar(_vec, _col); }

  private:
    std::vector<double, aligned<double, OLKMC_EIGEN_ALIGN>> _vec{};
    std::vector<Colour> _col{};
};

namespace periodic {

enum : std::size_t {
    x = 0b0001,
    y = 0b0010,
    z = 0b0100,
    w = 0b1000,

    non = 0,
};

}  // namespace periodic

// Provides details of the simulation supercell, all queries of the space in which the atoms exist
// are provided by this class. It is assumed all non-periodic atoms are within the Simbox extents.
class Simbox {
  public:
    Vec3<double> extents;
    Vec3<double> inv_extents;

    std::size_t periodic;

    Simbox() = default;

    Simbox(Vec3<double> const &extents, std::size_t periodic);

    // Maps atom into canonical cell, 0 <= r_i < extent_i for all i which are periodic. Non-periodic
    // atoms are within the simbox extents so v[i] * inv_extents less than 1 and  v[i] remains
    // unaffected, hence no non-periodic switch
    Vec3<double> canonicle_image(Vec3<double> const &v) const {
        return v - extents * (v * inv_extents).floor();
    }

    // Returns the smallest vector connecting a to a periodic image of b
    Vec3<double> min_image(Vec3<double> const &a, Vec3<double> const &b) const;

    template <class Archive> void serialize(Archive &ar) { ar(extents, inv_extents, periodic); }

    static Simbox load(toml::v2::table const &config);

  private:
    friend bool operator==(Simbox const &a, Simbox const &b) {
        return (a.extents == b.extents).all() && a.periodic == b.periodic;
    }
};

// Represents all the atoms in the simulation. Atoms are partitioned into "active" atoms which can
// move during the simulation and "boundary" atoms which remain fixed.
class Supercell : public Simbox {
  public:
    AtomVector activ;
    AtomVector bound;

    explicit Supercell(Simbox const &box) : Simbox{box} {}

    Supercell() = default;

    std::size_t size() const { return activ.size() + bound.size(); }

    // Compute the displacement of active from *this to other with the translational component
    // removed.
    VecN<double> active_disp(VecN<double> const &others) const;

    template <class Archive> void serialize(Archive &ar) {
        ar(cereal::base_class<Simbox>(this), activ, bound);
    }
};

// ----------------------------------------- //

// Augmentation of a Supercell that packages a supercell and a location to do SPS work
struct Workcell : Supercell {
    std::size_t centre;
};

// ----------------------------------------- //

// From LAMMPS compatible .xyz file
std::pair<Supercell, std::string> load_supercell(
    toml::v2::table const &config,
    std::unordered_map<std::string, std::uint16_t> const &);

// In LAMMPS compatible .xyz file
void dump_supercell(Supercell const &cell, std::string const &out_file, bool append = false);

// Automatic naming
void dump_supercell(Supercell const &cell);
