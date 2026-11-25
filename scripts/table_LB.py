"""LB(lam) at every eigenvalue lam <= X, in parallel; compare with data/table5000.txt (C++).
Usage: python scripts/table_LB.py 5000 [workers]"""
import sys, os
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
from squareecp import deficiency, lattice as L


def row(lam):
    lb = deficiency.lower_bound(lam)
    return lam, L.dim(lam), L.mult(lam), lb, lb - (L.mult(lam) - 1)


if __name__ == "__main__":
    X = int(sys.argv[1]); W = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    with Pool(W) as p:
        for r in p.imap(row, L.eigenvalues(X), chunksize=4):
            print(*r, flush=True)
