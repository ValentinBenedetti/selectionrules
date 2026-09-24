"""Nodes (closed submodels) of the N=2 minimal model k, their commuting lines and global indices.

Reproduces the numbers in N2_TQFT_note.pdf.  Conventions as in N2minimal.m / N2Gab.pdf:
fields (l,m,s), l+m+s even, (l,m,s) ~ (k-l, m+k+2, s+2), S-matrix eq. (2.6) with normalisation 1/(k+2).
Usage:  python3 n2_nodes.py 3 5 7 9 11
"""
import sys
from fractions import Fraction
import numpy as np


def key(f, k):
    K = k + 2
    redm = lambda m: ((m + K - 1) % (2 * K)) - (K - 1)
    reds = lambda s: ((s + 1) % 4) - 1
    l, m, s = f
    return min((l, redm(m), reds(s)), (k - l, redm(m + K), reds(s + 2)))


def build(k):
    K = k + 2
    seen, rs = set(), []
    for l in range(k + 1):
        for s in range(-1, 3):
            for m in range(-(k + 1), k + 3):
                if (l + m + s) % 2 or key((l, m, s), k) in seen:
                    continue
                seen.add(key((l, m, s), k))
                rs.append((l, m, s))
    rs.remove((0, 0, 0))
    rs = [(0, 0, 0)] + rs
    S = np.array([[np.sin(np.pi * (a[0] + 1) * (b[0] + 1) / K) * np.exp(1j * np.pi * a[1] * b[1] / K)
                   * np.exp(-1j * np.pi * a[2] * b[2] / 2) / K for b in rs] for a in rs])
    Nf = np.round(np.einsum('iu,ju,ku,u->ijk', S, S, S.conj(), 1 / S[0]).real).astype(int)
    return rs, S, Nf


def closure(gen, Nf):
    cur = set(gen) | {0}
    while True:
        new = set(cur)
        for a in cur:
            for b in cur:
                new |= set(np.nonzero(Nf[a, b])[0])
        if new == cur:
            return frozenset(cur)
        cur = new


def node(k, gens, data=None):
    """Smallest node containing gens; returns (fields, commuting lines, global index mu)."""
    rs, S, Nf = data or build(k)
    d = (S[0] / S[0, 0]).real
    T = closure([next(i for i, g in enumerate(rs) if key(g, k) == key(f, k)) for f in gens], Nf)
    F = [j for j in range(len(rs)) if all(abs(S[j, i] / S[0, i] - S[j, 0] / S[0, 0]) < 1e-8 for i in T)]
    mu = sum(d[j] ** 2 for j in F) ** 2
    return [rs[i] for i in sorted(T)], [rs[j] for j in F], mu


def h(f, k):
    l, m, s = f
    return (Fraction(l * (l + 2) - m * m, 4 * (k + 2)) + Fraction(s * s, 8)) % 1


if __name__ == "__main__":
    for k in map(int, sys.argv[1:] or [3, 5, 7, 9, 11]):
        K = k + 2
        dso3 = K / (4 * np.sin(np.pi / K) ** 2)
        uv = build(k)
        T, F, mu = node(k, [(0, -2, 0)], uv)
        T2, F2, mu2 = node(k, [(0, -2, 0), (k - 2, k - 2, 2)], uv)
        print(f"k={k}: dim so(3)_k={dso3:.4f}")
        print(f"  N_UV=<(0,-2,0)>            |T|={len(T):3d} |F|={len(F):3d} mu={mu:.4f}")
        print(f"  E_k=<(0,-2,0),(k-2,k-2,2)>  |T|={len(T2):3d} F={F2} mu={mu2:.4f}")
        print(f"  h(0,K,1) = {h((0, K, 1), k)}")
        if k - 2 >= 4:
            T3, F3, mu3 = node(k - 2, [(4, 2, 0)])
            print(f"  IR k'={k-2}: <(4,2,0)>        |T|={len(T3):3d} F={F3} mu={mu3:.4f}"
                  f"  h(0,K',1) = {h((0, K - 2, 1), k - 2)}")
