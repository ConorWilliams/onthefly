#include "local/environment.hpp"

#include <algorithm>

#include "config.hpp"
#include "local/geometry.hpp"
#include "utility.hpp"

namespace options {

Mechanism Mechanism::load(toml::v2::table const& config) {
    Mechanism opt;

    opt.r_tol = fetch<double>(config, "mechanism", "r_tol");
    opt.energy_abs_tol = fetch<double>(config, "mechanism", "abs_tol");
    opt.energy_frac_tol = fetch<double>(config, "mechanism", "frac_tol");
    opt.rel_cap_tol = fetch<double>(config, "mechanism", "rel_cap_tol");

    return opt;
}

}  // namespace options

bool MechBase::within_tol(MechBase const& other, double abs_tol, double frac_tol) const {
    double dif_activ = std::abs(other.activ_energy - activ_energy);

    if (dif_activ / activ_energy > frac_tol && dif_activ > abs_tol) {
        return false;
    }

    double dif_delta = std::abs(other.delta_energy - delta_energy);

    if (dif_delta / std::abs(delta_energy) > frac_tol && dif_delta > abs_tol) {
        return false;
    }

    return true;
}

bool Mechanism::within_tol(Mechanism const& other,
                           double abs_tol,
                           double frac_tol,
                           double r_tol) const {
    if (MechBase::within_tol(other, abs_tol, frac_tol)) {
        double sum_sq = 0;

        for (std::size_t i = 0; i < disp.size(); i++) {
            sum_sq += norm_sq(disp[i] - other.disp[i]);
        }

        return sum_sq < r_tol * r_tol;
    } else {
        return false;
    }
}

bool ProtoMech::within_tol(ProtoMech const& other,
                           double abs_tol,
                           double frac_tol,
                           double r_tol) const {
    return MechBase::within_tol(other, abs_tol, frac_tol) && norm(other.disp - disp) < r_tol;
}

std::size_t ProtoMech::find_centre() const {
    // Sanity
    CHECK(disp.size() >= 3, "");
    CHECK(disp.size() % 3 == 0, "");

    std::size_t centre = 0;
    double max = 0;

    for (int i = 0; i < disp.size() / 3; ++i) {
        double norm_sq = disp[3 * i + 0] * disp[3 * i + 0] + disp[3 * i + 1] * disp[3 * i + 1]
                         + disp[3 * i + 2] * disp[3 * i + 2];

        if (norm_sq > max) {
            max = norm_sq;
            centre = i;
        }
    }

    return centre;
}

bool Environment::try_push_mech(Mechanism&& m, double abs_tol, double frac_tol, double r_tol) {
    auto it = std::find_if(mechs.begin(), mechs.end(), [&](Mechanism const& other) {
        return m.within_tol(other, abs_tol, frac_tol, r_tol);
    });

    if (it == mechs.end()) {
        // std::cout << "New mechanism active=" << m.activ_energy << '\n';
        mechs.push_back(std::move(m));
        return true;
    } else {
        // std::cout << "Repeated mechanism activ=" << m.activ_energy << '\n';
        return false;
    }
}
