#include "potentials/ADP/potential.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

#include "config.hpp"
#include "potentials/ADP/data.hpp"
#include "potentials/neigh_reduce.hpp"
#include "supercell.hpp"

PotentialADP::PotentialADP(TabADP &&data) : _data(std::make_shared<TabADP>(std::move(data))) {}

std::unique_ptr<PotentialBase> PotentialADP::clone() const {
    return std::make_unique<PotentialADP>(*this);
}

// Compute energy
double PotentialADP::energy(Supercell const &cell) {
    _reduce.load(rcut(), cell);

    double v_sum = 0;
    double f_sum = 0;
    double u_sum = 0;
    double w_sum = 0;
    double epsilon_sum = 0;

    for (auto a = _reduce.begin_activ(); a != _reduce.begin_bound(); ++a) {
        double rho = 0;
        Vec3<double> mu_dipole = Vec3<double>::Zero();
        Mat3<double> lambda_quadropole = Mat3<double>::Zero();

        _reduce.neigh_reduce(a, [&](auto b, double r, Vec3<double> const &dr) {
            v_sum += _data->v(a->col.atomic, b->col.atomic)(r);
            rho += _data->phi(b->col.atomic)(r);

            for (std::size_t i = 0; i < 3; i++) {
                mu_dipole(i) += _data->u(a->col.atomic, b->col.atomic)(r) * (-dr[i]);
                for (std::size_t j = 0; j < 3; j++) {
                    lambda_quadropole(i, j)
                        += _data->w(a->col.atomic, b->col.atomic)(r) * dr[i] * dr[j];
                }
            }
        });

        f_sum += _data->f(a->col.atomic)(rho);
        u_sum += norm_sq(mu_dipole);
        w_sum += lambda_quadropole.squaredNorm();
        epsilon_sum += lambda_quadropole.trace() * lambda_quadropole.trace();
    }

    return (0.5 * v_sum) + f_sum + (0.5 * u_sum) + (0.5 * w_sum) - (epsilon_sum / 6);
}

// Compute force
void PotentialADP::gradient(Supercell const &cell, VecN<double> &out) {
    _reduce.load(rcut(), cell);

    // First sum computes density, dipole density u and quadrupole density v at each atom, runs over
    // active+boundary atoms
    for (auto a = _reduce.begin_activ(); a != _reduce.begin_ghost(); ++a) {
        double rho = 0;
        Vec3<double> mu_dipole = Vec3<double>::Zero();
        Mat3<double> lambda_quadropole = Mat3<double>::Zero();

        // Computes rho at atom
        _reduce.neigh_reduce(a, [&](auto b, double r, Vec3<double> const &dr) {
            rho += _data->phi(b->col.atomic)(r);

            for (std::size_t i = 0; i < 3; i++) {
                mu_dipole(i) += _data->u(a->col.atomic, b->col.atomic)(r) * (-dr[i]);
                for (std::size_t j = 0; j < 3; j++) {
                    lambda_quadropole(i, j)
                        += _data->w(a->col.atomic, b->col.atomic)(r) * dr[i] * dr[j];
                }
            }
        });

        // Compute F'(rho) at atom
        a->density = _data->f(a->col.atomic).grad(rho);
        a->dipole_density = mu_dipole;  // Note in NiCr ADP, this is exactly 0 as u(r) = 0
        a->quadropole_density = lambda_quadropole;
    }

    // Update density of ghost atoms
    _reduce.broadcast_ghost_data();

    out.resize(3 * cell.activ.size());  // Usually a no-op

    // Counter must be declared outside loop
    std::size_t k = 0;
    // Second sum computes force, only runs over active atoms
    for (auto i = _reduce.begin_activ(); i != _reduce.begin_bound(); ++i) {
        out[3 * k + 0] = 0;
        out[3 * k + 1] = 0;
        out[3 * k + 2] = 0;

        // Finds dr = r^{ji}
        _reduce.neigh_reduce(i, [&](auto j, double r, Vec3<double> const &dr) {
            // EAM terms
            double mag = _data->v(i->col.atomic, j->col.atomic).grad(r)
                         + i->density * _data->phi(j->col.atomic).grad(r)
                         + j->density * _data->phi(i->col.atomic).grad(r);
            mag /= r;
            for (std::size_t g = 0; g < 3; g++) {
                out[3 * k + g] -= mag * (-dr[g]);
                // Term 1, Acta Materialia 53 (2005) 4029–4041. Note = 0 in NiCr ADP as u(r) = 0
                out[3 * k + g] -= (i->dipole_density(g) - j->dipole_density(g))
                                  * _data->u(i->col.atomic, j->col.atomic)(r);
                // Term 5
                out[3 * k + g] += (1.0 / 3.0) * (-dr[g])
                                  * (i->quadropole_density.trace() + j->quadropole_density.trace())
                                  * (_data->w(i->col.atomic, j->col.atomic).grad(r) * r
                                     + 2 * _data->w(i->col.atomic, j->col.atomic)(r));
                for (std::size_t a = 0; a < 3; a++) {
                    // Term 2. Note = 0 in NiCr ADP as u(r) = 0
                    out[3 * k + g] -= (dr[a] * dr[g] / r)
                                      * (i->dipole_density(a) - j->dipole_density(a))
                                      * _data->u(i->col.atomic, j->col.atomic).grad(r);
                    // Term 3
                    out[3 * k + g] -= 2 * (-dr[a])
                                      * (i->quadropole_density(a, g) + j->quadropole_density(a, g))
                                      * _data->w(i->col.atomic, j->col.atomic)(r);
                    for (std::size_t b = 0; b < 3; b++) {
                        // Term 4
                        out[3 * k + g]
                            -= (-dr[a] * dr[b] * dr[g] / r)
                               * (i->quadropole_density(a, b) + j->quadropole_density(a, b))
                               * _data->w(i->col.atomic, j->col.atomic).grad(r);
                    }
                }
            }
        });

        ++k;
    }
}
