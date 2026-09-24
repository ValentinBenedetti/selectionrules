"""N=2-closed nodes of the N=2 minimal model k (odd k) for the flow k -> k-2.

A node is N=2 closed iff it contains the diagonal field (0,0,2) = G x Gbar, i.e. with every (l,m,s)
it contains its superpartner sector (l,m,s+2).  Reproduces the numbers of N2_TQFT_note.pdf.
Usage:  python3 n2_susy.py 3 5 7 9 11
"""
import sys
from fractions import Fraction
import numpy as np
from n2_nodes import build, closure, key


def index_of(f, k, rs):
    return next(i for i, g in enumerate(rs) if key(g, k) == key(f, k))


def analyze(k, gens, data):
    rs, S, Nf = data
    d = (S[0] / S[0, 0]).real
    T = closure([index_of(f, k, rs) for f in gens], Nf)
    F = [j for j in range(len(rs)) if all(abs(S[j, i] / S[0, i] - S[j, 0] / S[0, 0]) < 1e-8 for i in T)]
    return T, F, sum(d[j] ** 2 for j in F) ** 2


def supermultiplets(k, T, rs):
    """Superconformal primaries [l,m] (|m|<=l, l=m mod 2) whose NS module (l,m,0)+(l,m,2) lies in T."""
    keys = {key(rs[i], k) for i in T}
    out = []
    for l in range(k + 1):
        for m in range(-l, l + 1, 2):
            if key((l, m, 0), k) in keys:
                assert key((l, m, 2), k) in keys
                out.append((l, m))
    return out


if __name__ == "__main__":
    for k in map(int, sys.argv[1:] or [3, 5, 7, 9, 11]):
        K = k + 2
        data = build(k)
        rs = data[0]
        NS = {i for i, f in enumerate(rs) if f[2] % 2 == 0}
        print(f"=== k={k}  (c={Fraction(3 * k, K)}),  #fields={len(rs)}, #NS fields={len(NS)}")
        # minimal N=2-closed node of the perturbation
        T, F, mu = analyze(k, [(0, -2, 0), (0, -2, 2)], data)
        print(f" N^S_UV=<(0,-2,0),(0,-2,2)>: |T|={len(T)} contains (0,0,2): {index_of((0,0,2),k,rs) in T}"
              f"  |F|={len(F)} F={[rs[j] for j in F]}  mu={mu:.4f}")
        print("   supermultiplets [l,m] (h, q=m/K, Delta_bottom=2h, Delta_(s=2)=2h+1):")
        for (l, m) in supermultiplets(k, T, rs):
            h = Fraction(l * (l + 2) - m * m, 4 * K)
            tag = []
            if 2 * h < 2 and (l, m) != (0, 0):
                tag.append("bottom relevant")
            if 2 * h + 1 < 2 and (l, m) != (0, 0):
                tag.append("F-term relevant (N=2 preserving)")
            print(f"     [{l},{m:+d}] h={str(h):>6} q={str(Fraction(m, K)):>6}  {', '.join(tag)}")
        # massless flow: add F-term of X^{k-2}
        T2, F2, mu2 = analyze(k, [(0, -2, 0), (0, -2, 2), (k - 2, k - 2, 2)], data)
        print(f" M^S_UV=<N^S,(k-2,k-2,2)>: |T|={len(T2)} equals NS sector: {set(T2) == NS}"
              f"  F={[rs[j] for j in F2]} mu={mu2:.4f}")
        # N=2-closed nodes containing the perturbation: F-subrings of F(N^S) (checked by closure)
        d = (data[1][0] / data[1][0, 0]).real
        so3 = [index_of((L, 0, 0), k, rs) for L in range(0, k + 1, 2)]
        gen_all = all(set(closure([j], data[2])) >= set(so3) for j in so3[1:])
        print(f" every nontrivial so(3)_k line generates all of so(3)_k: {gen_all}")
        # IR
        if k - 2 >= 4:
            kp = k - 2
            dI = build(kp)
            T3, F3, mu3 = analyze(kp, [(4, 2, 0), (0, 0, 2)], dI)
            NSI = {i for i, f in enumerate(dI[0]) if f[2] % 2 == 0}
            print(f" IR k'={kp}: <(4,2,0),(0,0,2)>: |T|={len(T3)} equals NS: {set(T3) == NSI}"
                  f"  F={[dI[0][j] for j in F3]} mu={mu3:.4f}")
