"""Nodal domains of linear combinations of Dirichlet eigenfunctions of the square.

Code for the paper "The extended Courant property, Part I: nodal domains of linear
combinations of eigenfunctions on the square".

Modules
-------
lattice      support S(lam), kappa, multiplicity, lattice-polygon primitives
bound        the bound U(lam) of Theorem A (exact, by exhaustive search)
deficiency   the lower bound LB(lam) <= dim E_lam - U(lam) (exact dynamic programme)
hull_units   H(lam), the hull-unit count of Theorem C
nodal        numerical nodal-domain counts for explicit polynomials (evidence only)
"""
__version__ = "1.0.0"
