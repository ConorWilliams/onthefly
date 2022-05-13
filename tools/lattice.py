import numpy as np
import itertools


class atom:
    def __init__(self, col, vec):
        self.col = col
        self.vec = np.asarray(vec)


bcc = [
    atom(0, [0.0, 0.0, 0.0]),
    atom(0, [0.5, 0.5, 0.5]),
]


a = np.asarray([1, 0, 0])
b = np.asarray([0, 1, 0])
c = np.asarray([0, 0, 1])


def build(shape, param=2.8557, basis=[a, b, c], motif=bcc, key=lambda vec:  True):
    lat = []

    r1 = [i for i in range(shape[0])]
    r2 = [i for i in range(shape[1])]
    r3 = [i for i in range(shape[2])]

    for (a, b, c) in itertools.product(r1, r2, r3):
        point = a*basis[0] + b*basis[1] + c*basis[2]
        for at in motif:
            vec = param*(point+at.vec)
            if(key(vec)):
                lat.append(atom(at.col, vec))

    return param * np.asarray(shape), lat


def vacancy_at(points, tol=0.5):
    def near_any(vec):
        for p in points:
            if np.linalg.norm(vec - np.asarray(p)) < tol:
                return False
        return True

    return near_any


temp = 300.0

V = 11.64012 + 9.37798e-5 * temp + 3.643134e-7 * temp**2 - \
    1.851593e-10 * temp**3 + 5.669148e-14 * temp**4

a = (2*V) ** (1/3)


extents, lat = build([6, 6, 6], param=a
                     # ,key=vacancy_at(
                     #     [[a, a, a], [a, 2*a, a], [2*a, a, a], [2*a, 2*a, a], [1.5*a, 1.5*a, 1.5*a]]

                     # )
                     )


for at in lat:
    at.vec += 1.0

with open("../data/FeH/V0H1.xyz", 'w') as file:
    print(len(lat), file=file)
    print("Generated .xyz @", temp, "kelvin, extents =", extents, file=file)
    for at in lat:
        print(at.col, *at.vec, file=file)
