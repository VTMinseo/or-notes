"""Shared toy integer program for the Integer Programming basics pages.

    maximize   2 x1 + 7 x2
    subject to -5 x1 + 5 x2 <= 14
                7 x1 + 6 x2 <= 36
                0 <= x1 <= 5,  0 <= x2 <= 6,   x1, x2 integer

LP relaxation optimum (1.477, 4.277), value 32.89.  Integer optimum (2, 3), value 25.
Rounding the LP optimum down gives (1, 4), which is infeasible.
"""
from __future__ import annotations

import itertools
import math

import numpy as np
from scipy.optimize import linprog

c = np.array([2.0, 7.0])
A = np.array([[-5.0, 5.0], [7.0, 6.0]])
b = np.array([14.0, 36.0])
SUB = "₁₂"
LO = np.array([0.0, 0.0])
HI = np.array([5.0, 6.0])


def all_rows(A=A, b=b, lo=LO, hi=HI):
    """A x <= b together with the bounds, as one system G x <= h."""
    G = np.vstack([A, -np.eye(2), np.eye(2)])
    h = np.concatenate([b, -np.asarray(lo, float), np.asarray(hi, float)])
    return G, h


def lp(lo=LO, hi=HI, A=A, b=b, c=c, extra=()):
    """LP relaxation over the box [lo, hi] (and extra cuts (g, beta)). Returns (x, z) or (None, None)."""
    AA, bb = A, b
    if extra:
        AA = np.vstack([A] + [np.atleast_2d(g) for g, _ in extra])
        bb = np.concatenate([b, [beta for _, beta in extra]])
    if np.any(np.asarray(lo) > np.asarray(hi) + 1e-9):
        return None, None
    r = linprog(-c, A_ub=AA, b_ub=bb, bounds=list(zip(lo, hi)), method="highs")
    return (r.x, -r.fun) if r.status == 0 else (None, None)


def integer_points(A=A, b=b, lo=LO, hi=HI):
    pts = [(i, j) for i in range(int(lo[0]), int(hi[0]) + 1)
           for j in range(int(lo[1]), int(hi[1]) + 1)
           if np.all(A @ np.array([i, j]) <= b + 1e-9)]
    return np.array(pts, float)


def ip_optimum(A=A, b=b, lo=LO, hi=HI, c=c):
    P = integer_points(A, b, lo, hi)
    k = int(np.argmax(P @ c))
    return P[k], float(P[k] @ c)


def polygon(G, h, box=50):
    """Vertices (counter-clockwise) of {x : G x <= h} in the plane."""
    from _lib.orplot import clip_halfplane
    P = [(-box, -box), (box, -box), (box, box), (-box, box)]
    for g, v in zip(G, h):
        P = clip_halfplane(P, g, v)
    return P


# ---------------------------------------------------------------- branch and bound
def branch_and_bound(rule="best"):
    """LP-based branch and bound, branching on the most fractional variable.

    rule: "best" (open node with the best LP bound first) or "depth" (newest first).
    Returns a list of node dicts in the order they were processed: id, parent, depth,
    lo, hi, x, z, status ('branched', 'integer', 'pruned: bound', 'infeasible'),
    branch (var, value), the incumbent value after processing, and a label for the
    branching constraint that created the node.
    """
    nodes, out = [], []
    open_ = [dict(id=0, parent=None, depth=0, lo=LO.copy(), hi=HI.copy(), how="root")]
    incumbent, inc_x, next_id = -math.inf, None, 1
    while open_:
        # evaluate LP for open nodes lazily
        for nd in open_:
            if "z" not in nd:
                nd["x"], nd["z"] = lp(nd["lo"], nd["hi"])
        if rule == "best":
            open_.sort(key=lambda nd: -(nd["z"] if nd["z"] is not None else -math.inf))
            nd = open_.pop(0)
        else:
            nd = open_.pop()
        x, z = nd["x"], nd["z"]
        rec = dict(id=nd["id"], parent=nd["parent"], depth=nd["depth"], lo=nd["lo"].tolist(),
                   hi=nd["hi"].tolist(), x=None if x is None else x.tolist(), z=z, how=nd["how"],
                   branch=None)
        if x is None:
            rec["status"] = "infeasible"
        elif z <= incumbent + 1e-9:
            rec["status"] = "pruned: bound"
        else:
            frac = np.abs(x - np.round(x))
            if frac.max() < 1e-7:
                rec["status"] = "integer"
                if z > incumbent:
                    incumbent, inc_x = z, np.round(x)
            else:
                j = int(np.argmax(frac))
                v = x[j]
                rec["status"] = "branched"
                rec["branch"] = (j, float(v))
                h1 = nd["hi"].copy(); h1[j] = math.floor(v)
                l2 = nd["lo"].copy(); l2[j] = math.ceil(v)
                kids = [dict(id=next_id, parent=nd["id"], depth=nd["depth"] + 1, lo=nd["lo"].copy(),
                             hi=h1, how=f"x{SUB[j]} ≤ {math.floor(v)}"),
                        dict(id=next_id + 1, parent=nd["id"], depth=nd["depth"] + 1, lo=l2,
                             hi=nd["hi"].copy(), how=f"x{SUB[j]} ≥ {math.ceil(v)}")]
                next_id += 2
                open_ += kids
        rec["incumbent"] = None if incumbent == -math.inf else incumbent
        rec["incumbent_x"] = None if inc_x is None else inc_x.tolist()
        rec["open_bound"] = max([o.get("z") or -math.inf for o in open_] + [-math.inf])
        out.append(rec)
    return out


# ---------------------------------------------------------------- Chvátal–Gomory cuts
def cg_cut(xstar, rows, grid=(0, 1, 2, 3, 4, 5, 6), denoms=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)):
    """Search small rational multipliers u >= 0 for the Chvátal–Gomory cut
    floor(u G) x <= floor(u h) most violated by xstar. rows = (G, h) including bounds."""
    G, h = rows
    best = None
    m = len(G)
    for idx in itertools.combinations(range(m), 2):
        for num in itertools.product(grid, repeat=2):
            for d in denoms:
                u = np.zeros(m)
                u[list(idx)] = np.array(num) / d
                if not u.any():
                    continue
                a = np.floor(u @ G + 1e-9)
                beta = math.floor(u @ h + 1e-9)
                if not a.any():
                    continue
                viol = (a @ xstar - beta) / np.linalg.norm(a)
                if viol > 1e-6 and (best is None or viol > best[0] + 1e-9):
                    best = (viol, a, float(beta), u.copy())
    return best


def cutting_plane_rounds(k=4):
    """Repeatedly add the most violated CG cut found by cg_cut. Returns list of
    (x_lp, z_lp, cut or None)."""
    G, h = all_rows()
    cuts, out = [], []
    for _ in range(k):
        x, z = lp(extra=cuts)
        best = cg_cut(x, (G, h))
        out.append(dict(x=x, z=z, cut=None if best is None else (best[1], best[2]),
                        u=None if best is None else best[3]))
        if best is None or np.allclose(x, np.round(x)):
            break
        cuts.append((best[1], best[2]))
        G = np.vstack([G, best[1]]); h = np.append(h, best[2])
    return out


# ---------------------------------------------------------------- Lagrangian relaxation
RELAX = 1   # relax the second constraint, 7 x1 + 6 x2 <= 36


def lagrangian_value(lam):
    """L(lam) = max over integer x in X of c x + lam (b_R - a_R x), where X keeps the
    other constraint and the bounds. Returns (value, maximizer)."""
    keep = [i for i in range(len(A)) if i != RELAX]
    X = integer_points(A[keep], b[keep])
    vals = X @ c + lam * (b[RELAX] - X @ A[RELAX])
    k = int(np.argmax(vals))
    return float(vals[k]), X[k]


def lagrangian_dual(lams=None):
    lams = np.linspace(0, 2.0, 801) if lams is None else lams
    vals = np.array([lagrangian_value(l)[0] for l in lams])
    k = int(np.argmin(vals))
    return lams, vals, float(lams[k]), float(vals[k])


# ---------------------------------------------------------------- facility location
def facility_location(disaggregated, integer=False, m=4, n=8, seed=3):
    """Uncapacitated facility location with m sites and n customers (fixed seed).

    aggregated:     sum_j x_ij <= n y_i      (one constraint per site)
    disaggregated:  x_ij <= y_i              (one constraint per site-customer pair)
    Returns (objective, y, positions of sites, positions of customers, x matrix)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    rng = np.random.default_rng(seed)
    F = rng.uniform(0, 10, (m, 2)); C = rng.uniform(0, 10, (n, 2))
    f = np.array([30, 25, 35, 28.0])[:m]
    cost = np.linalg.norm(F[:, None] - C[None], axis=2) * 3
    N = m + m * n
    obj = np.concatenate([f, cost.ravel()])
    rows, hi = [], []
    for j in range(n):                       # each customer served once
        r = np.zeros(N); r[m + np.arange(m) * n + j] = 1; rows.append(r); hi.append(1)
    ub = np.array(hi, float); lb = ub.copy()
    link = []
    if not disaggregated:
        for i in range(m):
            r = np.zeros(N); r[m + i * n:m + (i + 1) * n] = 1; r[i] = -n; link.append(r)
    else:
        for i in range(m):
            for j in range(n):
                r = np.zeros(N); r[m + i * n + j] = 1; r[i] = -1; link.append(r)
    A_ = np.vstack(rows + link)
    lb = np.concatenate([lb, -np.inf * np.ones(len(link))])
    ub = np.concatenate([ub, np.zeros(len(link))])
    integ = np.r_[np.ones(m), np.zeros(m * n)] if integer else np.zeros(N)
    r = milp(obj, constraints=LinearConstraint(A_, lb, ub), bounds=Bounds(0, 1), integrality=integ)
    return float(r.fun), r.x[:m], F, C, r.x[m:].reshape(m, n)


def tree_layout(nodes):
    """x positions for a branch-and-bound tree (leaves spaced evenly, parents centred)."""
    kids = {}
    for nd in nodes:
        if nd["parent"] is not None:
            kids.setdefault(nd["parent"], []).append(nd["id"])
    for k in kids:
        kids[k].sort()
    pos, nxt = {}, [0]

    def place(v):
        if v not in kids:
            pos[v] = nxt[0]; nxt[0] += 1
        else:
            for w in kids[v]:
                place(w)
            pos[v] = sum(pos[w] for w in kids[v]) / len(kids[v])
    place(0)
    return pos


def bb_snapshots(rule="best"):
    """For each step k of branch_and_bound: the node processed, the set of open nodes
    afterwards, the incumbent and the global upper bound max(incumbent, open LP bounds)."""
    nodes = branch_and_bound(rule)
    byid = {nd["id"]: nd for nd in nodes}
    done, snaps = [], []
    for k, nd in enumerate(nodes):
        done.append(nd["id"])
        created = {m["id"] for m in nodes if m["parent"] in done}
        open_ = sorted(created - set(done))
        bounds = [byid[i]["z"] for i in open_ if byid[i]["z"] is not None]
        inc = nd["incumbent"]
        ub = max(bounds + ([inc] if inc is not None else []), default=None)
        snaps.append(dict(step=k, node=nd["id"], done=list(done), open=open_,
                          incumbent=inc, upper=ub))
    return nodes, snaps


def lagrangian_set():
    """X: integer points of the box satisfying the constraints that are kept."""
    keep = [i for i in range(len(A)) if i != RELAX]
    return integer_points(A[keep], b[keep])


def lagrangian_pieces():
    """Each x in X gives the line lam -> c x + lam (b_R - a_R x). Returns (X, intercepts, slopes)."""
    X = lagrangian_set()
    return X, X @ c, b[RELAX] - X @ A[RELAX]


def subgradient(iters=30, lam0=0.0, step0=0.08):
    """Projected subgradient method for min_{lam >= 0} L(lam) with steps step0 / sqrt(k + 1).
    Returns list of (lam, L(lam), maximizer x, subgradient g)."""
    lam, out = lam0, []
    for k in range(iters):
        val, x = lagrangian_value(lam)
        g = float(b[RELAX] - A[RELAX] @ x)
        out.append((lam, val, x.tolist(), g))
        lam = max(0.0, lam - step0 / np.sqrt(k + 1) * g)
    return out
