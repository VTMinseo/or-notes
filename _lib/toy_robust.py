"""Shared toy models for the Robust Optimization pages.

Everything follows the budget-of-uncertainty model of Bertsimas and Sim (2004),
"The Price of Robustness", Operations Research 52(1):35-53.

* A 2-variable LP with two uncertain constraints (figures and the interactive page).
* The bound B(n, Gamma) of their Theorem 3 and the simpler exp bound of Theorem 2.
* A random 0-1 knapsack instance in the spirit of their Section 6.1, solved for a
  range of Gamma with scipy's MILP solver (HiGHS).
"""
from __future__ import annotations

from math import comb, exp, floor

import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

from _lib.orplot import clip_halfplane

# ---------------------------------------------------------------------------
# 2-D toy:  max 3 x1 + 4 x2
#   row 1:  (2 ± 0.6) x1 + (1 ± 0.5) x2 <= 10
#   row 2:  (1 ± 0.4) x1 + (3 ± 0.9) x2 <= 15
#   x >= 0.  Both coefficients of each row are uncertain, so |J_i| = 2.
# ---------------------------------------------------------------------------
C = np.array([3.0, 4.0])
ABAR = np.array([[2.0, 1.0], [1.0, 3.0]])
AHAT = np.array([[0.6, 0.5], [0.4, 0.9]])
B = np.array([10.0, 15.0])
BOX = [(0.0, 0.0), (8.0, 0.0), (8.0, 8.0), (0.0, 8.0)]   # drawing window


def budget_vertices(n: int, gamma: float):
    """Vertices of {z in [0,1]^n : sum z <= gamma} (n = 2 here, enough for figures)."""
    assert n == 2
    g = min(max(gamma, 0.0), 2.0)
    if g <= 1:
        return [np.array([0.0, 0.0]), np.array([g, 0.0]), np.array([0.0, g])]
    return [np.array([0.0, 0.0]), np.array([1.0, 0.0]), np.array([0.0, 1.0]),
            np.array([1.0, g - 1.0]), np.array([g - 1.0, 1.0])]


def protection(x, ahat, gamma):
    """beta(x, Gamma): the |_Gamma_| largest a_hat_j |x_j| plus a fraction of the next."""
    v = np.sort(np.abs(np.asarray(ahat) * np.asarray(x)))[::-1]
    k = int(floor(gamma))
    return float(v[:k].sum() + (gamma - k) * (v[k] if k < len(v) else 0.0))


def row_halfplanes(i: int, gamma: float):
    """Linear pieces of  abar.x + beta(x, Gamma) <= b  on x >= 0 (one per vertex z)."""
    return [(ABAR[i] + AHAT[i] * z, B[i]) for z in budget_vertices(2, gamma)]


def robust_region(gammas):
    """Robust feasible polygon for per-row budgets (Gamma_1, Gamma_2)."""
    V = BOX
    for i, g in enumerate(gammas):
        for a, b in row_halfplanes(i, g):
            V = clip_halfplane(V, a, b)
    out = []
    for v in V:                                  # drop repeated vertices
        if not out or np.hypot(v[0] - out[-1][0], v[1] - out[-1][1]) > 1e-9:
            out.append(v)
    if len(out) > 1 and np.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1]) <= 1e-9:
        out.pop()
    return out


def uncertainty_set(i: int, gamma: float):
    """U_i(Gamma) = {abar + ahat * eta : |eta_j| <= 1, sum |eta_j| <= Gamma}, as a polygon."""
    g = min(max(gamma, 0.0), 2.0)
    if g == 0:
        etas = [(0.0, 0.0)]
    elif g <= 1:
        etas = [(g, 0), (0, g), (-g, 0), (0, -g)]
    else:
        h = g - 1
        etas = [(1, h), (h, 1), (-h, 1), (-1, h), (-1, -h), (-h, -1), (h, -1), (1, -h)]
    return [tuple(ABAR[i] + AHAT[i] * np.array(e, dtype=float)) for e in etas]


def robust_lp(gammas):
    """Solve the robust counterpart of Theorem 1 of Bertsimas-Sim as one LP.

    Variables: x (2), y (2), z_i (2), p_ij (2x2).  Returns (x, objective).
    """
    n, m = 2, 2
    nv = n + n + m + m * n
    ix = lambda j: j
    iy = lambda j: n + j
    iz = lambda i: 2 * n + i
    ip = lambda i, j: 2 * n + m + i * n + j
    A, b = [], []
    for i in range(m):
        row = np.zeros(nv)
        row[[ix(0), ix(1)]] = ABAR[i]
        row[iz(i)] = gammas[i]
        for j in range(n):
            row[ip(i, j)] = 1.0
        A.append(row); b.append(B[i])
        for j in range(n):                      # z_i + p_ij >= ahat_ij y_j
            row = np.zeros(nv)
            row[iz(i)] = -1; row[ip(i, j)] = -1; row[iy(j)] = AHAT[i, j]
            A.append(row); b.append(0.0)
    for j in range(n):                          # -y_j <= x_j <= y_j
        row = np.zeros(nv); row[ix(j)] = 1; row[iy(j)] = -1; A.append(row); b.append(0.0)
        row = np.zeros(nv); row[ix(j)] = -1; row[iy(j)] = -1; A.append(row); b.append(0.0)
    c = np.zeros(nv); c[:n] = -C
    res = linprog(c, A_ub=np.array(A), b_ub=np.array(b), bounds=[(0, None)] * nv,
                  method="highs")
    return res.x[:n], -res.fun


# ---------------------------------------------------------------------------
# Probability bounds (Section 4 of Bertsimas-Sim)
# ---------------------------------------------------------------------------
def bound_exp(n: int, gamma: float) -> float:
    """Theorem 2: exp(-Gamma^2 / (2n))."""
    return exp(-gamma ** 2 / (2 * n))


def bound_B(n: int, gamma: float) -> float:
    """Theorem 3: B(n, Gamma), tight for symmetric +-1 deviations."""
    nu = (gamma + n) / 2
    f = int(floor(nu))
    mu = nu - f
    s1 = sum(comb(n, l) for l in range(f, n + 1))
    s2 = sum(comb(n, l) for l in range(f + 1, n + 1))
    return ((1 - mu) * s1 + mu * s2) / 2 ** n


def gamma_for(n: int, eps: float) -> float:
    """Smallest Gamma with B(n, Gamma) <= eps (bisection; B decreases in Gamma)."""
    if bound_B(n, n) > eps:
        return float(n)
    lo, hi = 0.0, float(n)
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if bound_B(n, mid) > eps else (lo, mid)
    return hi


# ---------------------------------------------------------------------------
# Robust 0-1 knapsack, in the spirit of Section 6.1 (our own random instance)
# ---------------------------------------------------------------------------
def knapsack_instance(N: int = 200, seed: int = 7):
    rng = np.random.default_rng(seed)
    w = rng.integers(20, 30, size=N).astype(float)     # nominal weights 20..29
    c = rng.integers(16, 78, size=N).astype(float)     # values 16..77
    what = 0.1 * w                                     # 10% deviation
    return c, w, what, 4000.0


def robust_knapsack(c, w, what, cap, gamma):
    """max c.x  s.t.  w.x + Gamma z + sum p <= cap,  z + p_j >= what_j x_j,  x binary."""
    N = len(c)
    nv = N + 1 + N                                    # x, z, p
    obj = np.concatenate([-c, [0.0], np.zeros(N)])
    rows, lo, hi = [], [], []
    r = np.concatenate([w, [gamma], np.ones(N)]); rows.append(r); lo.append(-np.inf); hi.append(cap)
    M = np.zeros((N, nv))
    M[np.arange(N), np.arange(N)] = what
    M[:, N] = -1.0
    M[np.arange(N), N + 1 + np.arange(N)] = -1.0
    A = np.vstack([np.array(rows), M])
    lo = np.concatenate([lo, -np.inf * np.ones(N)])
    hi = np.concatenate([hi, np.zeros(N)])
    integrality = np.concatenate([np.ones(N), np.zeros(N + 1)])
    ub = np.concatenate([np.ones(N), np.inf * np.ones(N + 1)])
    res = milp(obj, constraints=LinearConstraint(A, lo, hi), integrality=integrality,
               bounds=Bounds(np.zeros(nv), ub), options={"mip_rel_gap": 1e-6})
    x = np.round(res.x[:N])
    return x, float(c @ x)


def violation_rate(x, w, what, cap, n_sims: int = 20000, kind: str = "pm1", seed: int = 0):
    """Monte-Carlo probability that the realized weight exceeds the capacity."""
    rng = np.random.default_rng(seed)
    idx = np.flatnonzero(x)
    if kind == "pm1":
        eta = rng.choice([-1.0, 1.0], size=(n_sims, len(idx)))
    else:
        eta = rng.uniform(-1.0, 1.0, size=(n_sims, len(idx)))
    load = w[idx].sum() + eta @ what[idx]
    return float(np.mean(load > cap + 1e-9))
