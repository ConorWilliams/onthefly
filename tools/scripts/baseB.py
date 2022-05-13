import math
import scipy.constants as const
import matplotlib.pyplot as plt

# heavyside step


Z_Fe = 26.0
Z_H = 1

elementary = 1.602176634e-19
e_0 = 8.854187812813e-12

r_bohr = const.physical_constants["Bohr radius"][0] * 1e10  # angstoms

r_s = (
    0.88534 * r_bohr / math.sqrt(pow(Z_H, 2 / 3) + pow(Z_Fe, 2 / 3))
)  # angstroms


# a_phi = [
#     14.0786236789212005,
#     -4.4526835887173704,
#     5.5025121262565992,
#     -1.0687489808214079,
#     -0.3461498208163201,
#     -0.0064991947759021,
#     -0.0357435602984102,
# ]

a_r_phi = [  # new
    (14.0786236766230779, 1.6),
    (-4.4526835638887965, 1.7),
    (5.5025349784052979, 1.8),
    (-1.0687331741292405, 2.0),
    (-0.3461226670484926, 2.5),
    (-0.0064991313802717, 3.2),
    (-0.0357322844877736, 4.2),
]

a_F = [  # new
    -0.0581047132616673,
    0.0022873205657864,
    -0.0000313966169286,
    0.0000013788174098,
    -0.0000000253074673,
    0.0000000001487789,
]


# r_phi = [1.6, 1.7, 1.8, 2.0, 2.5, 3.2, 4.2]


B = [  # new
    1242.154614241987,
    -6013.4610429013765,
    12339.275191444543,
    -12959.339514470237,
    6817.662603221567,
    -1422.130403271231,
]


def H(x):
    if x <= 0:
        return 0
    else:
        return 1


def screen(x):
    return (
        0.1818 * math.exp(-3.2 * x)
        + 0.5099 * math.exp(-0.9423 * x)
        + 0.2802 * math.exp(-0.4029 * x)
        + 0.02817 * math.exp(-0.2016 * x)
    )


r_s2 = (
    0.88534 * r_bohr / math.sqrt(pow(Z_Fe, 2 / 3) + pow(Z_Fe, 2 / 3))
)  # angstroms

B2 = [7.4122709384068, -0.64180690713367, -2.6043547961722, 0.6262539393123]


a_phi_Fe = [  # new
    (-27.444805994228, 2.2),
    (15.738054058489, 2.3),
    (2.2077118733936, 2.4),
    (-2.4989799053251, 2.5),
    (4.2099676494795, 2.6),
    (-0.77361294129713, 2.7),
    (0.80656414937789, 2.8),
    (-2.3194358924605, 3.0),
    (2.6577406128280, 3.3),
    (-1.0260416933564, 3.7),
    (0.35018615891957, 4.2),
    (-0.058531821042271, 4.7),
    (-0.0030458824556234, 5.3),
]


def phi_Fe_Fe(r):
    if r < 1.0:  # good
        return 9734.2365892908 * screen(r / r_s2)
    elif r <= 2.05:  # good
        return r * math.exp(sum(b * pow(r, i) for i, b in enumerate(B2)))
    else:  # good
        return r * sum(a * pow(r_p - r, 3) * H(r_p - r) for a, r_p in a_phi_Fe)


def phi_Fe_H(r):  # new, good
    if r < 0.6:
        return (
            screen(r / r_s)
            * Z_Fe
            * Z_H
            * const.e
            * 1e10
            / (4 * const.pi * const.epsilon_0)
        )
    elif r <= 1.2:
        return r * sum(b * pow(r, i) for i, b in enumerate(B))
    else:
        return r * sum(a * pow(r_p - r, 3) * H(r_p - r) for a, r_p in a_r_phi)


def F_Fe_rho(p):  # new good
    return (
        -math.sqrt(p)
        - 6.7314115586063e-4 * pow(p, 2)
        + 7.6514905604792e-8 * pow(p, 4)
    )


def F_H_rho(p):  # new, good
    return sum(a * pow(p, i) for i, a in enumerate(a_F, 1))


C_PHH = 1800  # new [same]


def rho_H_H_r(r):  # new, [near] good
    return C_PHH * pow(r, 2) * math.exp(-2 * r / r_bohr) * f_cut(r)


def s(r):
    return 0.5 * (1 - math.tanh(25 * (r - 0.9)))


r_cut_phi_HH = 2.4  # new


def f_cut(r):
    if r < r_cut_phi_HH:
        return math.exp(1 / (r - r_cut_phi_HH))
    else:
        return 0


r_0 = 0.74  # angstrom
lambd = 0.4899


def a_t(r):
    return (r - r_0) / (r_0 * lambd)


E_b = 2.37  # ev/atom


def E_mol(r):
    return -2 * E_b * (1 + a_t(r)) * math.exp(-a_t(r))


C_1_phi_HH = 0.0
C_2_phi_HH = 0.0


def phi_H_H(r):  # new {as good as rho_H_H}
    if r < r_cut_phi_HH:
        return r * (s(r) * (E_mol(r) - 2 * F_H_rho(rho_H_H_r(r))))
    else:
        return 0


def F_rho(i, rho):  # meta
    if i == "H":
        return F_H_rho(rho)
    elif i == "Fe":
        return F_Fe_rho(rho)
    else:
        return 40


a_p_Fe_H = [  # new
    (10.0073629218346891, 1.6),
    (32.4862873850836635, 1.8),
    (-0.9494211670931015, 2.0),
    (11.6683860903729624, 2.4),
    (-0.0147079871493827, 3.2),
    (0.4945807618408609, 4.2),
]


a_p_H_Fe = [  # new
    (11.1667357634216433, 1.5),
    (-3.0351469477486712, 2.0),
    (3.6092404272928578, 2.5),
    (0.0212508491354509, 3.0),
    (0.0303904795842773, 4.2),
]


a_p_Fe_Fe = [  # new
    (11.686859407970, 2.4),
    (-0.01471074009883, 3.2),
    (0.47193527075943, 4.2),
]


def rho(i, j, r):
    if i == "Fe" and j == "Fe":  # new, good
        return sum(a * pow(r_p - r, 3) * H(r_p - r) for a, r_p in a_p_Fe_Fe)
    if i == "Fe" and j == "H":  # new, good
        return sum(a * pow(r_p - r, 3) * H(r_p - r) for a, r_p in a_p_Fe_H)
    if i == "H" and j == "Fe":  # new, good
        return sum(a * pow(r_p - r, 3) * H(r_p - r) for a, r_p in a_p_H_Fe)
    if i == "H" and j == "H":  # meta
        return rho_H_H_r(r)
    else:
        print(i, j)
        return 41


def phi(i, j, r):  # meta
    if (i == "Fe" and j == "H") or (i == "H" and j == "Fe"):
        return phi_Fe_H(r)
    if i == "H" and j == "H":
        return phi_H_H(r)
    if i == "Fe" and j == "Fe":
        return phi_Fe_Fe(r)
    else:
        print(i, j)
        return 42


# x = [i * (6 / 1000) for i in range(1, 1000)]
# y = [rho("H", "H", r) / r for r in x]
#
#
# plt.plot(x, y)
#
# plt.show()
