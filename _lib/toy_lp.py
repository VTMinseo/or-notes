"""Shared toy LP for the Linear Programming pages, so every page shows the same numbers.

A workshop makes chairs (x1) and tables (x2):

    maximize   3 x1 + 4 x2            (profit)
    subject to   x1 + 2 x2 <= 8       (wood)
               2 x1 +   x2 <= 7       (labor hours)
               x1, x2 >= 0

Optimum x* = (2, 3), profit 18. Shadow prices y* = (5/3, 2/3).
The feasible region is the same polygon as on the Foundations page
"Polyhedra and Their Vertices".
"""
from __future__ import annotations

import itertools

import numpy as np
from scipy.optimize import linprog

A = np.array([[1.0, 2.0], [2.0, 1.0]])
b = np.array([8.0, 7.0])
c = np.array([3.0, 4.0])
RESOURCES = ["wood", "labor"]
PRODUCTS = ["chairs", "tables"]


def all_constraints(A=A, b=b):
    """Rows of A x <= b plus the nonnegativity rows -x <= 0."""
    n = A.shape[1]
    return np.vstack([A, -np.eye(n)]), np.concatenate([b, np.zeros(n)])


def vertices(A=A, b=b, tol=1e-9):
    """Vertices of {x >= 0 : A x <= b} in the plane, counter-clockwise."""
    G, h = all_constraints(A, b)
    pts = []
    for i, j in itertools.combinations(range(len(G)), 2):
        M = G[[i, j]]
        if abs(np.linalg.det(M)) < tol:
            continue
        v = np.linalg.solve(M, h[[i, j]])
        if np.all(G @ v <= h + 1e-7):
            pts.append(v)
    if not pts:
        return np.zeros((0, 2))
    pts = np.unique(np.round(pts, 9), axis=0)
    ctr = pts.mean(axis=0)
    ang = np.arctan2(pts[:, 1] - ctr[1], pts[:, 0] - ctr[0])
    return pts[np.argsort(ang)]


def solve(c=c, A=A, b=b):
    """Solve max c.x s.t. A x <= b, x >= 0 with HiGHS. Returns x, z, y (duals >= 0)."""
    res = linprog(-np.asarray(c), A_ub=A, b_ub=b, bounds=[(0, None)] * A.shape[1],
                  method="highs")
    if res.status != 0:
        return None
    return dict(x=res.x, z=-res.fun, y=-res.ineqlin.marginals)


def simplex_path(c=c, A=A, b=b):
    """Textbook tableau simplex (largest-coefficient rule) from the origin.

    Returns one dict per visited vertex: point, objective, basis names, reduced costs
    of the nonbasic variables, and the entering/leaving variables of the next pivot.
    """
    m, n = A.shape
    names = [f"x{j + 1}" for j in range(n)] + [f"s{i + 1}" for i in range(m)]
    T = np.hstack([A, np.eye(m), b[:, None]]).astype(float)
    obj = np.concatenate([-np.asarray(c, float), np.zeros(m + 1)])  # row: z - c.x = 0
    basis = list(range(n, n + m))
    steps = []
    while True:
        x = np.zeros(n + m)
        x[basis] = T[:, -1]
        red = {names[j]: float(-obj[j]) for j in range(n + m) if j not in basis}
        step = dict(point=x[:n].copy(), z=float(obj[-1]), basis=[names[k] for k in basis],
                    reduced=red, enter=None, leave=None)
        steps.append(step)
        j = int(np.argmin(obj[:-1]))
        if obj[j] >= -1e-12:
            break
        col = T[:, j]
        ratios = np.where(col > 1e-12, T[:, -1] / np.where(col > 1e-12, col, 1), np.inf)
        r = int(np.argmin(ratios))
        step["enter"], step["leave"] = names[j], names[basis[r]]
        T[r] /= T[r, j]
        for i in range(m):
            if i != r:
                T[i] -= T[i, j] * T[r]
        obj -= obj[j] * T[r]
        basis[r] = j
    return steps


def central_path(mus, c=c, A=A, b=b):
    """Points maximizing c.x + mu * sum(log slacks) (all constraints incl. x >= 0)."""
    G, h = all_constraints(A, b)
    x = np.array([1.0, 1.0])
    out = []
    for mu in mus:
        for _ in range(100):
            s = h - G @ x
            grad = c / mu - G.T @ (1 / s)
            H = -(G.T * (1 / s ** 2)) @ G
            dx = np.linalg.solve(H, -grad)
            t = 1.0
            while np.any(h - G @ (x + t * dx) <= 0):
                t *= 0.5
            x = x + t * dx
            if np.linalg.norm(t * dx) < 1e-10:
                break
        out.append(x.copy())
    return np.array(out)


def ellipsoid_iterates(G, h, center, radius, iters=8):
    """Central-cut ellipsoid method for finding a point of {x : G x <= h} in the plane.

    Returns a list of (center, Q) with E = {x : (x - center)^T Q^{-1} (x - center) <= 1};
    stops early once the center is feasible.
    """
    n = 2
    z = np.asarray(center, float)
    Q = radius ** 2 * np.eye(n)
    out = [(z.copy(), Q.copy())]
    for _ in range(iters):
        viol = G @ z - h
        if np.all(viol <= 1e-12):
            break
        a = G[int(np.argmax(viol))]
        g = Q @ a / np.sqrt(a @ Q @ a)
        z = z - g / (n + 1)
        Q = n ** 2 / (n ** 2 - 1) * (Q - 2 / (n + 1) * np.outer(g, g))
        out.append((z.copy(), Q.copy()))
    return out
