"""Recompute U(lam), kappa, N and LB at every eigenvalue lam in [lo, hi] (Table 1 of the paper).
Usage: python scripts/table_U.py 5 130 > data/U_table.txt      (lam = 130 takes several minutes)"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
from squareecp import bound, deficiency, hull_units, lattice as L

lo, hi = map(int, sys.argv[1:3])
print("lam r N kappa U N-U LB ceilH/2")
for lam in L.eigenvalues(hi):
    if lam < lo:
        continue
    U = bound.U(lam)
    N = L.dim(lam)
    print(lam, L.mult(lam), N, L.kappa(lam), U, N - U, deficiency.lower_bound(lam),
          (hull_units.H(lam) + 1) // 2, flush=True)
