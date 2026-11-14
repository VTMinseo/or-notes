"""A two-dimensional toy problem for the Dantzig–Wolfe pages.

    min  c^T x   s.t.  a^T x <= beta   (one complicating / linking constraint)
                       x in X          (an "easy" polytope)

with  c = (-2, -3),  a = (1, 2),  beta = 8,  and  X = conv(V)  for the vertices V
below. The vertices are deliberately not integer, so that X ∩ Z^2 has a smaller
convex hull (used for the integer-programming discussion).

Everything on the page (iterates, duals, bounds) is computed from this file.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linprog

from _lib.orplot import convex_hull, lattice_points

V = [(0.0, 0.0), (5.5, 0.0), (6.3, 2.2), (3.6, 4.7), (0.8, 3.6)]
C = np.array([-2.0, -3.0])
A_LINK = np.array([1.0, 2.0])
B_LINK = 8.0


def hrep(vertices):
    """{x : A x <= b} for a convex polygon given by its vertices."""
    H = convex_hull(vertices)
    A, b = [], []
    for i in range(len(H)):
        p, q = H[i], H[(i + 1) % len(H)]
        n = np.array([q[1] - p[1], p[0] - q[0]])
        A.append(n); b.append(n @ p)
    return np.array(A), np.array(b)


A_X, B_X = hrep(V)
INT_PTS = lattice_points(V, (0, 7), (0, 6))
A_XI, B_XI = hrep([tuple(p) for p in INT_PTS])


def lp_opt(A, b):
    """min c^T x over {A x <= b} ∩ {a^T x <= beta}."""
    res = linprog(C, A_ub=np.vstack([A, A_LINK]), b_ub=np.r_[b, B_LINK],
                  bounds=[(None, None)] * 2)
    return res.x, res.fun


def lp_relaxation():
    return lp_opt(A_X, B_X)


def dw_integer_bound():
    """Dantzig–Wolfe bound when pricing is over the integer points of X."""
    return lp_opt(A_XI, B_XI)


def integer_optimum():
    feas = [p for p in INT_PTS if A_LINK @ p <= B_LINK + 1e-9]
    best = min(feas, key=lambda p: C @ p)
    return np.array(best, float), float(C @ best)


def column_generation(first=(0.0, 0.0), max_iter=30):
    """Column generation for the Dantzig–Wolfe master of the toy problem.

    Restricted master:  min sum_v (c^T v) lam_v
                        s.t. sum_v (a^T v) lam_v <= beta   (dual pi <= 0)
                             sum_v lam_v = 1               (dual mu)
                             lam >= 0
    Pricing: min_{x in X} (c - pi a)^T x - mu  over the vertices of X.
    """
    cols = [np.array(first, float)]
    log = []
    for k in range(max_iter):
        Cm = np.array(cols)
        res = linprog(Cm @ C, A_ub=[Cm @ A_LINK], b_ub=[B_LINK],
                      A_eq=[np.ones(len(Cm))], b_eq=[1.0],
                      bounds=[(0, None)] * len(Cm), method="highs")
        pi = float(res.ineqlin.marginals[0]); mu = float(res.eqlin.marginals[0])
        x = Cm.T @ res.x
        pr = linprog(C - pi * A_LINK, A_ub=A_X, b_ub=B_X, bounds=[(None, None)] * 2)
        v = np.round(pr.x, 9) + 0.0
        rc = float(pr.fun - mu)
        lb = float(res.fun + min(rc, 0.0))
        stop = rc > -1e-9
        log.append(dict(k=k, cols=[c.copy() for c in cols], lam=res.x.copy(), x=x,
                        z=float(res.fun), pi=pi, mu=mu, new=None if stop else v,
                        rc=rc, lb=lb))
        if stop:
            break
        cols.append(v)
    return log
