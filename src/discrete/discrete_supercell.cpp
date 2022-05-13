#include "discrete/discrete_supercell.hpp"

#include "supercell.hpp"
#include "utility.hpp"

namespace options {

Discrete Discrete::load(toml::v2::table const &config) {
    Discrete opt;

    opt.r_neigh = config["discrete"]["r_neigh"].value_or(opt.r_neigh);
    opt.r_tol = config["discrete"]["r_tol"].value_or(opt.r_tol);

    opt.state_tol = config["discrete"]["state_tol"].value_or(opt.state_tol);
    opt.max_barrier = config["discrete"]["max_barrier"].value_or(opt.max_barrier);

    opt.lattice_type = config["discrete"]["lattice_type"].value_or(opt.lattice_type);
    opt.lattice_parameter = config["discrete"]["lattice_parameter"].value_or(opt.lattice_parameter);
    opt.repeat_units = config["discrete"]["repeat_units"].value_or(opt.repeat_units);

    opt.format = config["discrete"]["format"].value_or(opt.format);
    opt.fname = config["discrete"]["fname"].value_or(opt.fname);
    opt.load_from_disk = config["discrete"]["load_from_disk"].value_or(opt.load_from_disk);

    return opt;
}

}  // namespace options

void DiscreteSupercell::define_lattice(std::string const &lattice_type,
                                       Vec3<double> const &lattice_parameter,
                                       Vec3<double> const &repeat_units,
                                       Vec3<double> const &offset) {
    if (lattice_type == "FCC") {
        for (std::size_t i = 0; i < repeat_units[0]; i++) {
            for (std::size_t j = 0; j < repeat_units[1]; j++) {
                for (std::size_t k = 0; k < repeat_units[2]; k++) {
                    _lattice.emplace_back(
                        canonicle_image((Vec3<double>() << (i)*lattice_parameter[0] + offset[0],
                                         (j)*lattice_parameter[1] + offset[1],
                                         (k)*lattice_parameter[2] + offset[2])
                                            .finished()));
                    _lattice.emplace_back(canonicle_image(
                        (Vec3<double>() << (i + 0.5) * lattice_parameter[0] + offset[0],
                         (j + 0.5) * lattice_parameter[1] + offset[1],
                         (k)*lattice_parameter[2] + offset[2])
                            .finished()));
                    _lattice.emplace_back(
                        canonicle_image((Vec3<double>() << (i)*lattice_parameter[0] + offset[0],
                                         (j + 0.5) * lattice_parameter[1] + offset[1],
                                         (k + 0.5) * lattice_parameter[2] + offset[2])
                                            .finished()));
                    _lattice.emplace_back(canonicle_image(
                        (Vec3<double>() << (i + 0.5) * lattice_parameter[0] + offset[0],
                         (j)*lattice_parameter[1] + offset[1],
                         (k + 0.5) * lattice_parameter[2] + offset[2])
                            .finished()));
                }
            }
        }
    }
}

DiscreteSupercell load_discrete_supercell(Supercell const &init, toml::v2::table const &config) {
    DiscreteSupercell cell(init, options::Discrete::load(config));

    cell.define_lattice(cell._opt.lattice_type,
                        Vec3<double>::Constant(cell._opt.lattice_parameter),
                        Vec3<double>::Constant(cell._opt.repeat_units),
                        init.activ[0].vec);

    auto vacant = cell._lattice;
    cell._map.resize(cell._lattice.size());

    // Copy active atoms in order
    for (std::size_t i = 0; i < init.activ.size(); i++) {
        auto it
            = std::find_if(vacant.begin(), vacant.end(), [&](Vec3<double> const &lattice_point) {
                  return norm(init.min_image(init.activ[i].vec, lattice_point)) < cell._opt.r_tol;
              });
        if (it != vacant.end()) {
            auto it2 = std::find_if(
                cell._lattice.begin(), cell._lattice.end(), [&](Vec3<double> const &lattice_point) {
                    return (lattice_point[0] == (*it)[0]) && (lattice_point[1] == (*it)[1])
                           && (lattice_point[2] == (*it)[2]);
                });
            cell._map[it2 - cell._lattice.begin()] = i;
            cell.activ.emplace_back(it2 - cell._lattice.begin(), init.activ[i].col);
            vacant.erase(it);
        } else {
            ALWAYS_CHECK(false, "Error discretising lattice");
        }
    }

    for (std::size_t i = 0; i < init.bound.size(); i++) {
        auto it
            = std::find_if(vacant.begin(), vacant.end(), [&](Vec3<double> const &lattice_point) {
                  return norm(init.min_image(init.bound[i].vec, lattice_point)) < cell._opt.r_tol;
              });
        if (it != vacant.end()) {
            auto it2 = std::find_if(
                cell._lattice.begin(), cell._lattice.end(), [&](Vec3<double> const &lattice_point) {
                    return (lattice_point[0] == (*it)[0]) && (lattice_point[1] == (*it)[1])
                           && (lattice_point[2] == (*it)[2]);
                });
            cell._map[it2 - cell._lattice.begin()] = i + init.activ.size();
            cell.bound.emplace_back(it2 - cell._lattice.begin(), init.activ[i].col);
            vacant.erase(it);
        } else {
            ALWAYS_CHECK(false, "Error discretising lattice");
        }
    }

    for (std::size_t i = 0; i < vacant.size(); i++) {
        auto it = std::find_if(
            cell._lattice.begin(), cell._lattice.end(), [&](Vec3<double> const &lattice_point) {
                return (lattice_point[0] == vacant[i][0]) && (lattice_point[1] == vacant[i][1])
                       && (lattice_point[2] == vacant[i][2]);
            });
        cell._map[it - cell._lattice.begin()] = i + init.activ.size() + init.bound.size();
        cell.vacant.emplace_back(it - cell._lattice.begin(), Colour{0, Colour::vacant});
    }
    std::cout << vacant.size() << " Vacancies detected \n";
    return cell;
}

Supercell return_continuous_supercell(DiscreteSupercell const &init, bool include_vacancies) {
    Supercell cell(init);
    cell.activ.clear();
    cell.bound.clear();

    for (std::size_t i = 0; i < init.activ.size(); i++) {
        cell.activ.emplace_back(init._lattice[init.activ[i].lattice_site], init.activ[i].col);
    }
    for (std::size_t i = 0; i < init.bound.size(); i++) {
        cell.bound.emplace_back(init._lattice[init.bound[i].lattice_site], init.bound[i].col);
    }
    if (include_vacancies) {
        for (std::size_t i = 0; i < init.vacant.size(); i++) {
            cell.bound.emplace_back(init._lattice[init.vacant[i].lattice_site], init.vacant[i].col);
        }
    }
    return cell;
}

Supercell return_continuous_supercell(DiscreteSupercell const &init,
                                      Supercell const &match,
                                      std::unique_ptr<MinimiserBase> minimiser,
                                      std::unique_ptr<PotentialBase> ff,
                                      toml::v2::table const &config) {
    Supercell cell(init);
    cell.activ.clear();
    cell.bound.clear();

    for (std::size_t i = 0; i < init.activ.size(); i++) {
        cell.activ.emplace_back(init._lattice[init.activ[i].lattice_site], init.activ[i].col);
    }
    for (std::size_t i = 0; i < init.bound.size(); i++) {
        cell.bound.emplace_back(init._lattice[init.bound[i].lattice_site], init.bound[i].col);
    }
    if (!minimiser->minimise(cell, ff)) {
        ALWAYS_CHECK(false, "Minimisation failed on conversion from discrete cell to continuous");
    }

    DiscreteSupercell check = load_discrete_supercell(cell, config);
    CHECK(check.size() == init.size() && check.activ.size() == init.activ.size()
              && check.occupied() == init.occupied(),
          "Error in converting from discrete cell to continuous");
    for (std::size_t i = 0; i < init.size(); i++) {
        CHECK(check[i].col == init[i].col && check[i].lattice_site == init[i].lattice_site,
              "Error in converting from discrete cell to continuous");
    }
    CHECK(norm_sq(cell.active_disp(match.activ.view())) < init._opt.state_tol * init._opt.state_tol,
          "Error in converting from discrete cell to continuous"
              + std::to_string(norm_sq(cell.active_disp(match.activ.view()))));

    return cell;
}

void dump_discrete_supercell(DiscreteSupercell const &cell,
                             std::string const &out_file,
                             bool append) {
    std::ofstream outfile;
    if (append) {
        outfile.open(out_file, std::ios_base::app);
    } else {
        outfile.open(out_file);
    }

    outfile << cell.size() << std::endl;
    outfile << "Lattice=\"" << cell.extents[0] << " 0.0 0.0 0.0 " << cell.extents[1]
            << " 0.0 0.0 0.0 " << cell.extents[2] << "\"";

    outfile << std::setprecision(15);

    for (std::size_t i = 0; i < cell.activ.size(); ++i) {
        outfile << '\n' << cell.activ[i].col.atomic;
        outfile << ' ' << cell._lattice[cell.activ[i].lattice_site][0];
        outfile << ' ' << cell._lattice[cell.activ[i].lattice_site][1];
        outfile << ' ' << cell._lattice[cell.activ[i].lattice_site][2];
    }

    for (std::size_t i = 0; i < cell.bound.size(); ++i) {
        outfile << '\n' << cell.bound[i].col.atomic;
        outfile << ' ' << cell._lattice[cell.bound[i].lattice_site][0];
        outfile << ' ' << cell._lattice[cell.bound[i].lattice_site][1];
        outfile << ' ' << cell._lattice[cell.bound[i].lattice_site][2];
    }

    for (std::size_t i = 0; i < cell.vacant.size(); ++i) {
        outfile << '\n' << cell.vacant[i].col.atomic + 99;
        outfile << ' ' << cell._lattice[cell.vacant[i].lattice_site][0];
        outfile << ' ' << cell._lattice[cell.vacant[i].lattice_site][1];
        outfile << ' ' << cell._lattice[cell.vacant[i].lattice_site][2];
    }
    outfile << "\n";
}
