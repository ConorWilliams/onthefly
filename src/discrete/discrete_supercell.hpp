#pragma once
#include <vector>

#include "minimise/minimiser_base.hpp"
#include "potentials/potential_base.hpp"
#include "supercell.hpp"
#include "utility.hpp"

namespace options {

struct Discrete {
    double r_neigh = 3;  // (Angstrom), Nearest neighbour seperation a/sqrt(2) < r_neigh < a for fcc
    double r_tol = 0.1;  // Maximum distance of an atom in a supercell from a lattice point to
                         // consider the atom to lie on that lattice point

    double state_tol = 0.1;  // Norm squared difference between initial supercell and minimised
                             // discrete supercell to consider equivalent
    double max_barrier = 4.0;

    std::string lattice_type = "FCC";  // Details of the arrangements of lattice points
    double lattice_parameter = 3.52;
    double repeat_units = 10;

    std::string format = "portable_binary";  // Catalogue containing mechanisms acting on geometries
                                             // resolved on to lattice points
    std::string fname = "olkmc_discrete.cat";
    bool load_from_disk = false;

    static Discrete load(toml::v2::table const &config);
};
}  // namespace options

// Class resembling an AtomVector. However, the position is now discretised onto a fixed set of
// lattice sites
class DiscreteAtomVector {
  public:
    struct DiscreteAtomView {
        std::size_t &lattice_site;
        Colour &col;
    };

    struct ConstDiscreteAtomView {
        std::size_t const &lattice_site;
        Colour const &col;
    };

    std::size_t size() const { return _col.size(); }

    DiscreteAtomView operator[](std::size_t i) {
        CHECK(i < size(), "Bad access");
        return {_pos[i], _col[i]};
    }

    ConstDiscreteAtomView operator[](std::size_t i) const {
        CHECK(i < size(), "Bad access");
        return {_pos[i], _col[i]};
    }

    std::vector<std::size_t> &view() { return _pos; }

    std::vector<std::size_t> const &view() const { return _pos; }

    void emplace_back(std::size_t const &pos, Colour const &col) {
        _pos.insert(_pos.end(), pos);
        _col.insert(_col.end(), col);
    }

    void clear() {
        _pos.clear();
        _col.clear();
    }

  private:
    std::vector<std::size_t> _pos{};
    std::vector<Colour> _col{};
};

// Class resembling a Supercell, which tracks the current lattice sites of all atoms and stores the
// fixed positions of lattice points to reconstruct positions Also tracks vacant sites, as lattice
// sites without an atom on
class DiscreteSupercell : public Simbox {
  public:
    DiscreteAtomVector activ;
    DiscreteAtomVector bound;
    DiscreteAtomVector vacant;

    std::size_t size() const { return activ.size() + bound.size() + vacant.size(); }
    std::size_t occupied() const { return activ.size() + bound.size(); }

    DiscreteAtomVector::DiscreteAtomView operator[](std::size_t i) {
        CHECK(i < size(), "Bad access");
        if (i < activ.size()) {
            return activ[i];
        } else if (i < activ.size() + bound.size()) {
            return bound[i - activ.size()];
        } else {
            return vacant[i - activ.size() - bound.size()];
        }
    }

    DiscreteAtomVector::ConstDiscreteAtomView operator[](std::size_t i) const {
        CHECK(i < size(), "Bad access");
        if (i < activ.size()) {
            return activ[i];
        } else if (i < activ.size() + bound.size()) {
            return bound[i - activ.size()];
        } else {
            return vacant[i - activ.size() - bound.size()];
        }
    }

    explicit DiscreteSupercell(Simbox const &box, options::Discrete opt) : Simbox{box}, _opt{opt} {}

    DiscreteSupercell() = default;

    // List of positions of lattice points
    std::vector<Vec3<double>> _lattice{};
    // _map gives the index of a given atom in AtomVector such cell[_map[i]].lattice_site = i
    std::vector<std::size_t> _map{};

    void define_lattice(std::string const &lattice_type,
                        Vec3<double> const &lattice_parameter,
                        Vec3<double> const &repeat_units,
                        Vec3<double> const &offset);

    friend DiscreteSupercell load_discrete_supercell(Supercell const &cell,
                                                     toml::v2::table const &config);
    friend Supercell return_continuous_supercell(DiscreteSupercell const &init,
                                                 bool include_vacancies);
    friend Supercell return_continuous_supercell(DiscreteSupercell const &init,
                                                 Supercell const &match,
                                                 std::unique_ptr<MinimiserBase> minimiser,
                                                 std::unique_ptr<PotentialBase> ff,
                                                 toml::v2::table const &config);

    options::Discrete _opt;
};

// Constructs the lattice from config and compares to supercell to match each atom with a lattice
// point
DiscreteSupercell load_discrete_supercell(Supercell const &cell, toml::v2::table const &config);
// Unrelaxed i.e. atoms on lattice sites
Supercell return_continuous_supercell(DiscreteSupercell const &init,
                                      bool include_vacancies = false);
// Relaxes Discrete supercell, checks that it is the same as supercell match within tolerance
Supercell return_continuous_supercell(DiscreteSupercell const &init,
                                      Supercell const &match,
                                      std::unique_ptr<MinimiserBase> minimiser,
                                      std::unique_ptr<PotentialBase> ff,
                                      toml::v2::table const &config);

void dump_discrete_supercell(DiscreteSupercell const &cell,
                             std::string const &out_file,
                             bool append = false);
