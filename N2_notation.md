# Notation for the primary fields of the N = 2 minimal models (N2Gab.pdf)

Reference: F. Ambrosino, M. R. Gaberdiel, Y. Nakayama, *N = 2 RG flows, Non-Invertible
Symmetries and Matrix Factorisations*, arXiv:2608.02717, Section 2.1 (eqs. 2.1–2.7).

## 1. The model

The k-th N = 2 minimal model (k = 1, 2, 3, ...) has central charge

    c = 3k / (k+2)                                              (2.1)

and is realised as the coset

    (N=2)_k = [ su(2)_k ⊕ so(2)_1 ] / u(1)                      (2.2)

where so(2)_1 describes two free fermions. The u(1) in the denominator is u(1)_{2(k+2)}.
In the Landau–Ginzburg description it corresponds to the superpotential W = X^{k+2}
(the paper also uses d = k+2).

## 2. Labels of the primary fields: (l, m, s)

The irreducible representations are those of the **bosonic subalgebra** of the N = 2
superconformal algebra. Each one is labelled by three integers:

| label | origin | range |
|-------|--------|-------|
| `l` | su(2)_k highest weight (twice the spin) | l = 0, 1, ..., k |
| `m` | u(1)_{2(k+2)} charge | defined mod 2(k+2) |
| `s` | so(2)_1 label (free fermions) | defined mod 4 |

Constraints:

* **Selection rule:** `l + m + s` must be even. Triples with l+m+s odd are not allowed.
* **Field identification:**

      (l, m, s) ~ (k - l, m + k + 2, s + 2)

  Two triples related like this describe **the same field**. The identification has no
  fixed points, because it always shifts s by 2.

* `s` even (s = 0, 2) gives the **Neveu–Schwarz (NS)** sector. `s` odd (s = ±1) gives the
  **Ramond (R)** sector. s and s+2 are the two halves of a full N = 2 supermodule, i.e. the
  bosonic and fermionic parts of the same N = 2 representation.

### Counting

Before any constraint there are (k+1) · 2(k+2) · 4 triples. The parity rule halves this,
and the field identification halves it again. That leaves

    number of primaries = 2 (k+1)(k+2)

distinct fields of the bosonic algebra. For example, k = 1 gives 12, k = 2 gives 24, and
k = 3 gives 40.

Square brackets `[l, m, s]` denote the equivalence class under field identification.
Defects (topological lines) are written `D[L, M, S]`. They use the same labels because the
Cardy construction associates one defect with each primary.

## 3. Conformal weight and U(1) charge

For the highest-weight state of (l, m, s):

    h(l, m, s) = [ l(l+2) - m^2 ] / [ 4(k+2) ] + s^2/8     mod 1     (2.3)
    q(l, m, s) = s/2 - m/(k+2)                             mod 2     (2.4)

The paper writes q "mod 1". Because m is defined mod 2(k+2) and s mod 4, the formula
itself is well defined mod 2, and it is invariant under the field identification.

The formula for h is exact only when (l, m, s) lies in the "standard range" (roughly
|m - s| ≤ l). Otherwise it is correct only mod 1. That is enough for the T-matrix. For
example, `(0,0,2)` gives 1/2, but the true lowest weight of that module is 3/2 (it is G_{-3/2}|0⟩).

Useful special fields:

* `(0,0,0)` is the vacuum, and `(0,0,2)` holds the supercurrents G^±.
* `(l, ±l, 0)` are the (anti)chiral primaries, with h = l/(2(k+2)) and |q| = 2h.
* `(0,-2,0) ~ (k,k,2)` is the sector of the least relevant perturbation used in the paper
  (eq. 2.9).

## 4. Spectrum

The paper uses the charge-conjugation modular invariant:

    H = ⊕_{[l,m,s]}  H_(l,m,s) ⊗ H̄_(l,-m,-s)                    (2.5)

## 5. Modular S-matrix

    S_(l,m,s),(l',m',s') = N_k · sin[ π(l+1)(l'+1)/(k+2) ] · exp[ iπ m m'/(k+2) ] · exp[ -iπ s s'/2 ]     (2.6)

It is the product of the su(2)_k, u(1) and so(2)_1 S-matrices. Under the field
identification it picks up the factor (-1)^{l'+m'+s'} = 1, so it is well defined on
equivalence classes.

**Normalisation.** The paper prints N_k = 1/√(2(k+2)). On the 2(k+1)(k+2) inequivalent
fields that matrix is **not unitary**. The correct coset normalisation is

    N_k = 1/(k+2)

(= 2 × √(2/(k+2)) × 1/√(2(k+2)) × 1/2: the three factor normalisations times the orbit
factor 2). I checked numerically for k = 1, ..., 6 that with N_k = 1/(k+2):

* S S† = 1,
* S² = C, the charge conjugation (l,m,s) → (l,-m,-s),
* (S T)³ = S², with T = diag exp[2πi (h - c/24)],
* the Verlinde formula gives non-negative integer fusion coefficients.

Ratios such as the defect eigenvalues do not depend on N_k:

    γ^{[L,M,S]}_{[l,m,s]} = S_(L,M,S),(l,m,s) / S_(0,0,0),(l,m,s)          (2.7)

## 6. Example: k = 1 (c = 1), 12 fields

This is the order produced by `fieldsN2[1]` in `N2minimal.m`:

| # | representative (l,m,s) | identified with | h mod 1 | q |
|---|------|------|------|------|
| 1 | (0,0,0) | (1,3,2) | 0 | 0 |
| 2 | (0,-1,-1) | (1,2,1) | 1/24 | -1/6 |
| 3 | (0,1,-1) | (1,-2,1) | 1/24 | -5/6 |
| 4 | (0,3,-1) | (1,0,1) | 3/8 | -3/2 |
| 5 | (0,-2,0) | (1,1,2) | 2/3 | 2/3 |
| 6 | (0,2,0) | (1,-1,2) | 2/3 | -2/3 |
| 7 | (0,-1,1) | (1,2,-1) | 1/24 | 5/6 |
| 8 | (0,1,1) | (1,-2,-1) | 1/24 | 1/6 |
| 9 | (0,3,1) | (1,0,-1) | 3/8 | -1/2 |
| 10 | (0,-2,2) | (1,1,0) | 1/6 | 5/3 |
| 11 | (0,0,2) | (1,3,0) | 1/2 | 1 |
| 12 | (0,2,2) | (1,-1,0) | 1/6 | 1/3 |

## 7. How the Mathematica code chooses representatives

`N2minimal.m` works like this:

1. It lists every triple with l = 0..k, s ∈ {-1,0,1,2} and m ∈ {-(k+1), ..., k+2}
   (one period each), keeping only those with l+m+s even.
2. For each triple it builds a canonical key: it reduces both the triple and its image
   under the identification to these fixed ranges and takes the smaller of the two.
3. `DeleteDuplicatesBy` keeps the first triple for each key. Since l is scanned from 0
   upward, the kept representative has l ≤ k/2 whenever possible.
4. It moves the vacuum (0,0,0) to position 1. The routines in `NonUnitary.nb`
   (`verlinde`, `dimensions`, `comuttinglines`) assume index 1 is the identity.

`smat[[i,j]]` is then the S-matrix element between `rs[[i]]` and `rs[[j]]`.
