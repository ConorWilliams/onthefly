import numpy as np
import argparse

from datetime import datetime

from baseBB import F_rho, rho, phi

NPTS_RHO = int(10000)  # number of rho points, multiple of 5
DELTA_RHO = 2.9999999999999999e-02

R_CUT = 5.3000000000000000e00  # angstroms

NPTS_DIST = NPTS_RHO  # number of r points, multiple of 5
DELTA_DIST = R_CUT / NPTS_DIST

parser = argparse.ArgumentParser(description="Make an setfl EAM file")

parser.add_argument(
    "-o",
    "--out",
    metavar="Outfile",
    type=str,
    action="store",
    help="The outfile to write to",
)

args = parser.parse_args()

with open(args.out, "w") as f:
    print("Fe-H", file=f)
    print("C. J. Williams", file=f)
    print("Rendered at:", datetime.now(), file=f)

    print("2 Fe H", file=f)
    print(
        f"{NPTS_RHO} {DELTA_RHO:.16e} {NPTS_DIST} {DELTA_DIST:.16e} {R_CUT:.16e}",
        file=f,
    )
    print(" 26 55.847 2.855300 BCC", file=f)

    for j in range(NPTS_RHO // 5):
        f0 = F_rho("Fe", (j * 5 + 0) * DELTA_RHO)
        f1 = F_rho("Fe", (j * 5 + 1) * DELTA_RHO)
        f2 = F_rho("Fe", (j * 5 + 2) * DELTA_RHO)
        f3 = F_rho("Fe", (j * 5 + 3) * DELTA_RHO)
        f4 = F_rho("Fe", (j * 5 + 4) * DELTA_RHO)
        print(f"{f0:.16e}  {f1:.16e}  {f2:.16e}  {f3:.16e}  {f4:.16e}", file=f)

    for atom in ["Fe", "H"]:
        for j in range(NPTS_DIST // 5):
            p0 = rho("Fe", atom, (j * 5 + 0) * DELTA_DIST)
            p1 = rho("Fe", atom, (j * 5 + 1) * DELTA_DIST)
            p2 = rho("Fe", atom, (j * 5 + 2) * DELTA_DIST)
            p3 = rho("Fe", atom, (j * 5 + 3) * DELTA_DIST)
            p4 = rho("Fe", atom, (j * 5 + 4) * DELTA_DIST)
            print(
                f"{p0:.16e}  {p1:.16e}  {p2:.16e}  {p3:.16e}  {p4:.16e}", file=f
            )

    print("1 1.008 1.8 BCC", file=f)

    for j in range(NPTS_RHO // 5):
        f0 = F_rho("H", (j * 5 + 0) * DELTA_RHO)
        f1 = F_rho("H", (j * 5 + 1) * DELTA_RHO)
        f2 = F_rho("H", (j * 5 + 2) * DELTA_RHO)
        f3 = F_rho("H", (j * 5 + 3) * DELTA_RHO)
        f4 = F_rho("H", (j * 5 + 4) * DELTA_RHO)
        print(f"{f0:.16e}  {f1:.16e}  {f2:.16e}  {f3:.16e}  {f4:.16e}", file=f)

    for atom in ["Fe", "H"]:
        for j in range(NPTS_DIST // 5):
            p0 = rho("H", atom, (j * 5 + 0) * DELTA_DIST)
            p1 = rho("H", atom, (j * 5 + 1) * DELTA_DIST)
            p2 = rho("H", atom, (j * 5 + 2) * DELTA_DIST)
            p3 = rho("H", atom, (j * 5 + 3) * DELTA_DIST)
            p4 = rho("H", atom, (j * 5 + 4) * DELTA_DIST)
            print(
                f"{p0:.16e}  {p1:.16e}  {p2:.16e}  {p3:.16e}  {p4:.16e}", file=f
            )

    for a, b in [("Fe", "Fe"), ("Fe", "H"), ("H", "H")]:
        for j in range(NPTS_DIST // 5):
            v0 = phi(a, b, (j * 5 + 0) * DELTA_DIST)
            v1 = phi(a, b, (j * 5 + 1) * DELTA_DIST)
            v2 = phi(a, b, (j * 5 + 2) * DELTA_DIST)
            v3 = phi(a, b, (j * 5 + 3) * DELTA_DIST)
            v4 = phi(a, b, (j * 5 + 4) * DELTA_DIST)
            print(
                f"{v0:.16e}  {v1:.16e}  {v2:.16e}  {v3:.16e}  {v4:.16e}", file=f
            )
