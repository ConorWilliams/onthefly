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


a_phi = [
    14.0786236789212005,
    -4.4526835887173704,
    5.5025121262565992,
    -1.0687489808214079,
    -0.3461498208163201,
    -0.0064991947759021,
    -0.0357435602984102,
]

a_F = [
    -0.0581256120818134,
    0.0022854552833736,
    -0.0000314202805805,
    0.0000013764132084,
    -0.0000000253707731,
    0.0000000001483685,
]

r_phi = [1.6, 1.7, 1.8, 2.0, 2.5, 3.2, 4.2]

B = [
    1242.1709168218642,
    -6013.566711223783,
    12339.540893927151,
    -12959.66163724488,
    6817.850021676971,
    -1422.1723964897117,
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

B2 = [14.996917289290, -20.533174190155, 14.002591780752, -3.6473736591143]

a_phi_Fe = [
    (195.92322853994, 2.1),
    (17.516698453315, 2.2),
    (1.4926525164290, 2.3),
    (6.4129476125197, 2.4),
    (-6.8157461860553, 2.5),
    (9.6582581963600, 2.6),
    (-5.3419002764419, 2.7),
    (1.7996558048346, 2.8),
    (-1.4788966636288, 3.0),
    (1.8530435283665, 3.3),
    (-0.64164344859316, 3.7),
    (0.24463630025168, 4.2),
    (-0.057721650527383, 4.7),
    (0.023358616514826, 5.3),
    (-0.0097064921265079, 6.0),
]


# paper 2 equ 2
def phi_Fe_Fe(r):
    if r < 0.9:  # good
        return (
            screen(r / r_s2)
            * Z_Fe
            * Z_Fe
            * const.e
            * 1e10
            / (4 * const.pi * const.epsilon_0)
        )
    elif r <= 1.95:  # good
        return r * math.exp(sum(b * pow(r, i) for i, b in enumerate(B2)))
    else:  # good
        return r * sum(a * pow(r_p - r, 3) * H(r_p - r) for a, r_p in a_phi_Fe)


def phi_Fe_H(r):
    if r < 0.6:  # good
        return (
            screen(r / r_s)
            * Z_Fe
            * Z_H
            * const.e
            * 1e10
            / (4 * const.pi * const.epsilon_0)
        )
    elif r <= 1.2:  # good
        return r * sum(b * pow(r, i) for i, b in enumerate(B))
    else:  # good
        return r * sum(
            a * pow(r_p - r, 3) * H(r_p - r) for a, r_p in zip(a_phi, r_phi)
        )


def F_Fe_rho(p):  # good
    return -math.sqrt(p) - 0.00034906178363530 * pow(p, 2)


def F_H_rho(p):  # good
    return sum(a * pow(p, i) for i, a in enumerate(a_F, 1))


C_PHH = 1800


def rho_H_H_r(r):  # good
    return C_PHH * pow(r, 2) * math.exp(-2 * r / r_bohr) * f_cut(r)


def s(r):
    return 0.5 * (1 - math.tanh(25 * (r - 0.9)))


r_cut_phi_HH = 2.3


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


def phi_H_H(r):  # good
    if r <= r_cut_phi_HH:
        return r * (
            s(r) * (E_mol(r) - 2 * F_H_rho(rho_H_H_r(r)))
            + (1 - s(r)) * (C_1_phi_HH * f_cut(r) + C_2_phi_HH * rho_H_H_r(r))
        )
    else:
        return 0


def F_rho(i, rho):  # done
    if i == "H":
        return F_H_rho(rho)
    elif i == "Fe":
        return F_Fe_rho(rho)
    else:
        return 40


a_p_Fe_H = [
    (10.0073629216300581, 1.6),
    (32.4861983261490295, 1.8),
    (-0.9494226032063788, 2.0),
    (11.6659812262450338, 2.4),
    (-0.0147080251458273, 3.2),
    (0.4943383753319843, 4.2),
]

a_p_H_Fe = [
    (11.1667357634216433, 1.5),
    (-3.0351307365078730, 2.0),
    (3.6096144794370653, 2.5),
    (0.0212509034775648, 3.0),
    (0.0303914939946250, 4.2),
]

a_p_Fe_Fe = [
    (11.686859407970, 2.4),
    (-0.014710740098830, 3.2),
    (0.47193527075943, 4.2),
]


def rho(i, j, r):
    if i == "Fe" and j == "Fe":
        return sum(a * pow(r_p - r, 3) * H(r_p - r) for a, r_p in a_p_Fe_Fe)
    if i == "Fe" and j == "H":
        return sum(a * pow(r_p - r, 3) * H(r_p - r) for a, r_p in a_p_Fe_H)
    if i == "H" and j == "Fe":
        return sum(a * pow(r_p - r, 3) * H(r_p - r) for a, r_p in a_p_H_Fe)
    if i == "H" and j == "H":
        return rho_H_H_r(r)
    else:
        print(i, j)
        return 41


def phi(i, j, r):
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
