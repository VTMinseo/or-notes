"""A one-dimensional toy problem for the Benders decomposition pages.

    min  2.5 x + f(x),   x in {0, 1, 2, 3, 4}

where the subproblem (an LP in y, solved once x is fixed) is

    f(x) = min { y :  y >= 9 - 4x,  y >= 5 - x,  y >= 0 }.

Written as  min{ d^T y : B y >= b - A x }  with  d = 1,  B = (1, 1, 1)^T,
b = (9, 5, 0)^T,  A = (4, 1, 0)^T.  The dual feasible region is
Pi = { pi >= 0 : pi_1 + pi_2 + pi_3 = 1 },  whose extreme points are the unit
vectors; each gives one Benders optimality cut  t >= pi^T (b - A x):

    e1 -> t >= 9 - 4x,    e2 -> t >= 5 - x,    e3 -> t >= 0.

Everything on the pages (iterates, bounds, cuts) is computed from this file.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linprog

C = 2.5                       # first-stage cost per unit of x
XL, XU = 0.0, 4.0             # bounds on x
PIECES = [(9.0, 4.0), (5.0, 1.0), (0.0, 0.0)]   # cut t >= a - g x  for (a, g)
X_INT = np.arange(int(XL), int(XU) + 1)


def f(x: float) -> float:
    """Value function of the subproblem (max over dual extreme points)."""
    return max(a - g * x for a, g in PIECES)


def oracle(x: float):
    """Classical Benders oracle: solve the subproblem dual at x.

    Returns f(x) and the optimality cut (a, g) from the optimal dual extreme point.
    """
    vals = [a - g * x for a, g in PIECES]
    k = int(np.argmax(vals))
    return vals[k], PIECES[k]


def master_lp(cuts, extra=()):
    """LP master: min C x + t  s.t.  t >= a - g x for each cut, (extra) rows, x in [XL, XU].

    `extra` holds additional cuts written as (a, g) meaning t >= a - g x.
    """
    rows = list(cuts) + list(extra)
    A = [[-g, -1.0] for a, g in rows]
    b = [-a for a, g in rows]
    res = linprog([C, 1.0], A_ub=A, b_ub=b, bounds=[(XL, XU), (None, None)])
    return res.x, res.fun


def classic_benders(max_iter: int = 20):
    """Benders on the LP relaxation, starting from the trivial cut t >= 0.

    Returns a list of iterations with the master point, bounds and the cut added.
    """
    cuts = [(0.0, 0.0)]
    log, ub = [], np.inf
    for k in range(max_iter):
        (x, t), lb = master_lp(cuts)
        fx, cut = oracle(x)
        ub = min(ub, C * x + fx)
        stop = fx <= t + 1e-9
        log.append(dict(k=k, x=x, t=t, lb=lb, fx=fx, ub=ub,
                        cuts=list(cuts), new_cut=None if stop else cut))
        if stop:
            break
        cuts.append(cut)
    return log


def integer_optimum():
    vals = [C * x + f(x) for x in X_INT]
    k = int(np.argmin(vals))
    return X_INT[k], vals[k]


def routine1(xhat: float, that: float, split=(1.0, 2.0), init=((0.0, 0.0),), max_iter=20):
    """Routine 1 of Fang, Sin & Byeon (disjunctiveOracle) on the toy problem.

    Disjunction  x <= split[0]  OR  x >= split[1].  Normalization: reverse polar with
    direction (0, 1), so the LP finds the smallest tau >= 0 with
    (xhat, that + tau) in conv(P_1 ∪ P_2), using only the cuts known so far.
    Each iteration rescales the two disjunct points and asks `oracle` whether they
    lie in L; violated points produce new Benders cuts that are added for every term.

    Returns a list of iterations: tau, cut (alpha0, alpha_x, alpha_t) read off the
    duals, the rescaled points checked, and the cuts known at that iteration.
    """
    known = list(init)
    terms = [("le", split[0]), ("ge", split[1])]
    log = []
    for rho in range(max_iter):
        snapshot = list(known)
        n = 7  # (wx, wt, w0) for r = 0, 1, then tau
        A_ub, b_ub = [], []
        for r, (kind, val) in enumerate(terms):
            o = 3 * r
            for a, g in known:                         # wt >= a w0 - g wx
                row = np.zeros(n); row[o + 1] = -1; row[o + 2] = a; row[o] = -g
                A_ub.append(row); b_ub.append(0.0)
            row = np.zeros(n); row[o] = -1; row[o + 2] = XL; A_ub.append(row); b_ub.append(0.0)
            row = np.zeros(n); row[o] = 1; row[o + 2] = -XU; A_ub.append(row); b_ub.append(0.0)
            row = np.zeros(n)
            if kind == "le":
                row[o] = 1; row[o + 2] = -val
            else:
                row[o] = -1; row[o + 2] = val
            A_ub.append(row); b_ub.append(0.0)
        A_eq = [np.r_[0, 0, 1, 0, 0, 1, 0], np.r_[1, 0, 0, 1, 0, 0, 0], np.r_[0, 1, 0, 0, 1, 0, -1]]
        b_eq = [1.0, xhat, that]
        cost = np.zeros(n); cost[6] = 1.0
        res = linprog(cost, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                      bounds=[(None, None), (None, None), (0, None)] * 2 + [(0, None)],
                      method="highs")
        w, tau = res.x, res.x[6]
        a0, ax, at = res.eqlin.marginals            # cut:  a0 + ax x + at t <= 0
        checked, new = [], False
        for r in range(2):
            w0 = w[3 * r + 2]
            if w0 > 1e-7:
                xp, tp = w[3 * r] / w0, w[3 * r + 1] / w0
                fx, cut = oracle(xp)
                violated = tp < fx - 1e-9
                checked.append(dict(r=r, w0=w0, x=xp, t=tp, fx=fx, violated=violated, cut=cut))
                if violated and cut not in known:
                    known.append(cut); new = True
        log.append(dict(rho=rho, tau=tau, alpha=(a0, ax, at), checked=checked,
                        known=snapshot))
        if not new:
            break
    return log
