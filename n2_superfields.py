"""Nodes built from full N=2 superfields (superprimary together with all its superdescendants).

NS superfields [l,m] = (l,m,0)+(l,m,2), l+m even, [l,m] ~ [k-l, m+k+2].  For odd k every superfield has a
unique representative with l even (then m even, defined mod 2K).  Fusion: su(2)_k on l, addition on m.
Lines commuting with the N=2 chiral algebra are labelled by the same superfields; [L,M] acts on [l,m] by
(S_{Ll}/S_{0l}) exp(i pi M m / K), and the quantum dimension of [l,m] is the su(2)_k one.
Usage:  python3 n2_superfields.py 3 5 7 9 11
"""
import sys
import numpy as np


def ring(k):
    K = k + 2
    fields = [(l, m) for l in range(0, k + 1, 2) for m in range(0, 2 * K, 2)]
    return K, fields


def canon(l, m, k):
    K = k + 2
    m %= 2 * K
    if l % 2:                      # use [l,m] ~ [k-l, m+K] to get l even
        l, m = k - l, (m + K) % (2 * K)
    return (l, m)


def fuse(a, b, k):
    (l1, m1), (l2, m2) = a, b
    return {canon(l, m1 + m2, k) for l in range(abs(l1 - l2), min(l1 + l2, 2 * k - l1 - l2) + 1, 2)}


def closure(gens, k):
    cur = {(0, 0)} | {canon(*g, k) for g in gens}
    while True:
        new = set(cur)
        for a in cur:
            for b in cur:
                new |= fuse(a, b, k)
        if new == cur:
            return cur
        cur = new


def dsu2(l, k):
    return np.sin(np.pi * (l + 1) / (k + 2)) / np.sin(np.pi / (k + 2))


def lines_commuting(T, k):
    K, fields = ring(k)
    act = lambda L, M, l, m: (np.sin(np.pi * (L + 1) * (l + 1) / K) / np.sin(np.pi * (l + 1) / K)) \
        * np.exp(1j * np.pi * M * m / K)
    return [(L, M) for (L, M) in fields
            if all(abs(act(L, M, l, m) - dsu2(L, k)) < 1e-9 for (l, m) in T)]


def mu(F, k):
    return sum(dsu2(L, k) ** 2 for (L, M) in F) ** 2


if __name__ == "__main__":
    for k in map(int, sys.argv[1:] or [3, 5, 7, 9, 11]):
        K, fields = ring(k)
        dso3 = sum(dsu2(L, k) ** 2 for L in range(0, k + 1, 2))
        print(f"=== k={k}: {len(fields)} NS superfields (= (k+1)(k+2)/2), dim so(3)_k = {dso3:.4f}")
        T = closure([(k, k)], k)                      # superfield of the perturbation, X^k
        F = lines_commuting(T, k)
        print(f"  <[k,k]> = {sorted(T)}  (Z_K: [0,2n])")
        print(f"     F = {F}   mu = {mu(F, k):.4f}   dim^2 = {dso3**2:.4f}")
        T2 = closure([(k, k), (k - 2, k - 2)], k)     # add X^{k-2}
        F2 = lines_commuting(T2, k)
        print(f"  <[k,k],[k-2,k-2]>: {len(T2)} superfields (all: {len(T2) == len(fields)}), F = {F2}, mu = {mu(F2, k):.4f}")
        # every superfield node containing [k,k]: product of a subring of so(3)_k with Z_K
        so3_gen = all(len(closure([(L, 0)], k)) == (k + 1) // 2 for L in range(2, k + 1, 2))
        print(f"  each [L,0], L>0 even, generates all of so(3)_k: {so3_gen}")
        if k - 2 >= 4:
            kp = k - 2
            T3 = closure([(4, 2)], kp)
            F3 = lines_commuting(T3, kp)
            print(f"  IR k'={kp}: <[4,2]> has {len(T3)} of {len(ring(kp)[1])} superfields, F = {F3}, mu = {mu(F3, kp):.4f}")
