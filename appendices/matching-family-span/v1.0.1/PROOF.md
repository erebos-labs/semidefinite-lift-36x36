# Exact full-span obstruction for the published matching-test family

**Follow-up version 1.0.1 · 12 October 2026**

**Result:** the matching tests in equation (26), over their admissible domains
and with both prescribed rotation branches, together with the first-block
constant tests, span the entire 36-dimensional real test space. A finite
witness consists of four constants and 32 matching tests. In the published
ordered coefficient basis, its determinant is exactly **−1/8000**.

This completes the requested bounded investigation. The established matrix
order **36**, auxiliary count **24**, and representation of the entire closed
cone, including singular X, are preserved. No smaller lift is claimed.

## 1. Authority and scope

The authoritative source is
[erebos-labs/semidefinite-lift-36x36 at
2d2c170206ab937bcf56ff3290066cc3b2dfdd99](https://github.com/erebos-labs/semidefinite-lift-36x36/tree/2d2c170206ab937bcf56ff3290066cc3b2dfdd99),
the published v1.0.2 identified in the request. The following files were
retrieved at that exact commit and read:

- `THEOREM_PROOF.md`, especially Sections 3, 5, 6 and 8;
- `COORDINATES.md`;
- `generator.py`;
- `certificates/rank_kernel.json`.

Their unchanged bytes, source URLs and SHA-256 hashes accompany this handoff.
The retrieved repository-root MIT license is preserved as `source/LICENSE`.
The baseline's source-derived theorem exposition and adapted mathematical
notes retain Apache-2.0 terms; MIT does not relicense them. The Apache text,
component scope and attribution are in `LICENSES/Apache-2.0.txt`,
`LICENSING.md`, `NOTICE` and `ATTRIBUTION.md`. The source files retain some `v1.0.0` headers;
the requested commit, rather than those retained strings or an earlier local
draft, determines the definitions used here. This report's version number
versions only this separate follow-up.

The existing coefficient-rank and common-kernel certificates were read as
context. Their calculations were **not rerun**. This investigation forms and
checks only a matrix of test-vector coefficients. It does not form the 47
coefficient matrices of the semidefinite pencil.

## 2. The precise family whose span is computed

Let w = exp(i theta). The published test space is

\[
\mathcal W=\left\{
(\operatorname{Re}a,\operatorname{Im}a,
 \operatorname{Re}b,\operatorname{Im}b,u):
a=A+Bw^{-2},\ b=C+D_Fw^2
\right\},
\]

where A, B, C, D_F belong to C^4 and u belongs to R^4. The subscript F
distinguishes the fourth complex Fourier coefficient from a symmetric matrix
used below. Coordinates are ordered exactly as in the pinned source:

\[
(\Re A,\Im A,\Re B,\Im B,
  \Re C,\Im C,\Re D_F,\Im D_F,u),
\tag{F1}
\]

with coordinates 0, 1, 2, 3 within each four-vector group. Uniqueness of the
coefficients of 1 and w^{-2} in a, and of 1 and w^2 in b, proves that this is
a real basis of dimension 36.

Define the constant-test family

\[
\mathcal C=\{(v,0,0,0,0):v\in\mathbb R^4\}.
\tag{F2}
\]

These are the tests used in Section 8 of the source to infer X positive
semidefinite. Four basis vectors of this *linear family* will be selected
below. The assertion is not that four individual scalar inequalities suffice
to test positive semidefiniteness of X.

For the matching family, choose any real u, any real symmetric 4 by 4 matrix
D_sym, and

\[
z=\tau(c,s,0)+\rho(0,0,1),\qquad
\tau>0,\quad c^2+s^2=1,\quad \rho\in\mathbb R.
\tag{F3}
\]

Thus (z_1,z_2) is nonzero. Put w_0=c+is and

\[
U_0=\frac1{\sqrt2}
\begin{pmatrix}
-\bar w_0&-w_0\\
s\bar w_0^2&s w_0^2\\
c\bar w_0^2&-c w_0^2
\end{pmatrix},\qquad
V_0=\frac1{\sqrt2}
\begin{pmatrix}0&0\\i&i\\1&1\end{pmatrix},\qquad
R_0=\tau U_0+\rho V_0.
\tag{F4}
\]

For u=(u_0,u'), set

\[
J_u=(u'\ \ -u_0 I_3),\quad
S=J_uD_{\rm sym}J_u^T,\quad N=SR_0,
\]

and, numbering rows of N from 1 to 3, define

\[
v=N_3-iN_2,\quad
g=v_++\overline{v_-},\quad h=v_+-\overline{v_-},
\quad \delta=|g|^2+|h|^2.
\tag{F5}
\]

If delta is positive, use exactly the published rotation

\[
\kappa=w_0^2w^{-2},\quad
a_0=\frac{|g|^2+|h|^2\kappa}{\delta},\quad
\overline b_0=\frac{gh(1-\kappa)}{\delta},\quad
M=\begin{pmatrix}a_0&b_0\\-\overline b_0&\overline a_0\end{pmatrix}.
\tag{F6}
\]

If delta is zero, define M=I_2 directly. There is no division and no limiting
choice in this branch. Set script-C=R_0 M. Realify its two columns in the
order (Re column 1, Im column 1, Re column 2, Im column 2), calling the four
real columns b_1,...,b_4. The matching test is

\[
q(u,D_{\rm sym},z)=
(-D_{\rm sym}J_u^Tb_1,-D_{\rm sym}J_u^Tb_2,
 -D_{\rm sym}J_u^Tb_3,-D_{\rm sym}J_u^Tb_4,u).
\tag{F7}
\]

This is equation (26) of the pinned proof. Let \(\mathcal M\) denote all functions
(F7) with parameters (F3), arbitrary real u, arbitrary symmetric D_sym,
and their prescribed branch (F6) or M=I_2. In particular D_sym need not be
positive semidefinite.

## 3. Explicit coefficient-generation rule

The following rule expands every selected test exactly over Q(sqrt(2)). It
also displays containment in W for the entire admissible family.

Write

\[
-D_{\rm sym}J_u^TR_0=(p_1,p_2).
\]

For delta positive set

\[
\gamma=|g|^2/\delta,\qquad
\lambda=|h|^2/\delta,\qquad
\mu=gh/\delta,\qquad t=w_0^2.
\]

Then the four complex coefficient vectors of (F7) are

\[
\begin{aligned}
A&=\gamma p_1-\mu p_2,&
B&=t(\lambda p_1+\mu p_2),\\
C&=\overline\mu p_1+\gamma p_2,&
D_F&=\bar t(-\overline\mu p_1+\lambda p_2).
\end{aligned}
\tag{F8}
\]

For delta zero the rule is instead

\[
A=p_1,\quad B=0,\quad C=p_2,\quad D_F=0.
\tag{F9}
\]

These identities follow by multiplying the two columns in (F6); in
particular, conjugation changes w^{-2} to w^2. Formula (F8) does not treat
complex realification as an interchangeable coordinate order.

For exact rational computation, let Rhat=sqrt(2) R_0 and obtain ghat, hhat
using (F5) with Rhat. Then ghat=sqrt(2) g, hhat=sqrt(2) h and
deltahat=2 delta. The three ratios gamma, lambda, mu are unchanged. Computing
(F8) with (-D_sym J_u^T Rhat) therefore gives sqrt(2) times the first 32
real coefficient coordinates. The final four coordinates u are unchanged.

For the finite parameters in the next section, Rhat, ghat and hhat are
Gaussian rational. Thus if T is the actual 36 by 36 real coefficient matrix
of the selected functions, then

\[
V=\operatorname{diag}(\sqrt2 I_{32},I_4)T
\quad\text{has rational entries},\qquad
\det V=2^{16}\det T.
\tag{F10}
\]

## 4. The complete finite admissible witness

All indices on u and on 4 by 4 matrices in this section are **zero based**.
Write e_i for a standard basis vector, E_ii=e_i e_i^T,
F_ij=e_i e_j^T+e_j e_i^T, I=I_4, and 1=(1,1,1,1)^T.

The first four selected tests, in column order, are

\[
q_1=(e_0/\sqrt2,0,0,0,0),\quad\ldots,\quad
q_4=(e_3/\sqrt2,0,0,0,0).
\tag{F11}
\]

Every other column is exactly (F7), with tau=1 and the parameters below.
The table's delta is the *unscaled* source delta in (F5).

| Column | w_0 | rho | u | D_sym | delta |
|---:|:---:|---:|:---|:---|---:|
| 5 | 1 | 0 | e_0 | 0 | 0 |
| 6 | 1 | 0 | e_0 | I | 2 |
| 7 | 1 | 0 | e_0 | E_11 | 0 |
| 8 | 1 | 0 | e_0 | E_33 | 2 |
| 9 | 1 | 0 | e_0 | F_01 | 0 |
| 10 | 1 | 0 | e_0 | F_12 | 2 |
| 11 | 1 | 0 | e_0 | F_13 | 2 |
| 12 | 1 | 0 | e_0 | F_23 | 2 |
| 13 | 1 | 0 | e_0 | I+F_01 | 2 |
| 14 | 1 | 0 | e_0 | I+F_03 | 2 |
| 15 | 1 | 0 | e_0 | I+F_13 | 4 |
| 16 | 1 | 0 | e_0 | I+F_23 | 4 |
| 17 | 1 | 0 | e_1 | 0 | 0 |
| 18 | 1 | 0 | e_2 | 0 | 0 |
| 19 | 1 | 0 | e_3 | 0 | 0 |
| 20 | 1 | 0 | e_3 | I+F_02 | 2 |
| 21 | 1 | 0 | e_0+e_1 | I+F_23 | 4 |
| 22 | 1 | 0 | e_0+e_2 | I+F_03 | 4 |
| 23 | 1 | 0 | e_0+e_3 | I+F_01 | 10 |
| 24 | 1 | 0 | e_0+e_3 | I+F_02 | 10 |
| 25 | 1 | 0 | e_2+e_3 | I+F_01 | 4 |
| 26 | 1 | 0 | e_2+e_3 | I+F_03 | 4 |
| 27 | 1 | 0 | 1 | F_03 | 12 |
| 28 | 1 | 0 | 1 | I+F_12 | 20 |
| 29 | 1 | 1 | e_0 | I | 10 |
| 30 | 1 | 1 | e_0 | F_02 | 0 |
| 31 | 1 | 1 | e_0 | F_12 | 2 |
| 32 | 1 | 1 | e_0 | F_23 | 2 |
| 33 | i | 0 | e_0 | I | 2 |
| 34 | i | 0 | e_0 | F_01 | 0 |
| 35 | i | 0 | e_0 | F_12 | 2 |
| 36 | i | 0 | e_0 | F_13 | 2 |

This is an explicitly admissible list. Each matrix is symmetric, and some
F_ij are indefinite, as the source permits. Every tau is positive. The only
z values used are (1,0,0), (1,0,1), and (0,1,0), so the nonzero first-two-
coordinates hypothesis holds without any continuity argument. There are
eight matching tests in the delta-zero branch and 24 in the positive branch.

As a sign-and-order check, column 6 has u=e_0, D_sym=I, z=(1,0,0).
It gives delta=2, gamma=0, lambda=1, mu=0, and

\[
A=C=0,\quad B=(0,-1,0,1)^T/\sqrt2,\quad
D_F=(0,-1,0,-1)^T/\sqrt2,\quad u=e_0.
\]

Its only nonzero scaled coordinates are therefore
\(V_{9,5}=-1\), \(V_{11,5}=1\), \(V_{25,5}=-1\),
\(V_{27,5}=-1\), and \(V_{32,5}=1\). Both indices here are zero based,
as in the machine files; column index 5 is column 6 of the table.

## 5. Exact rank certificate

`independent/span_witness.json` stores all 1,296 rational entries of V, both
as rows and columns, and every parameter in the table as an exact rational
string. These entries are determined by (F4)–(F11); the saved matrix is not
an independent assumption.

Let d_j be the least common multiple of the denominators in column j of V.
In the table's one-based column order, the values different from 1 are

\[
\begin{aligned}
d_{16}=d_{21}=d_{22}=d_{25}=d_{26}&=2,\\
d_{23}=d_{24}=d_{28}=d_{29}&=5,\\
d_{27}&=3.
\end{aligned}
\tag{F12}
\]

Every other d_j is 1, so their product is \(2^5 5^4 3=60000\). Set

\[
W=V\operatorname{diag}(d_1,\ldots,d_{36})\in\mathbb Z^{36\times36}.
\]

`exact_rank_certificate.json` supplies this integer W and an integer matrix
U for which direct integer multiplication gives

\[
WU=UW=60I_{36}.
\tag{F13}
\]

The stored U has 330 nonzero entries, with absolute value at most 360.
Equation (F13) alone proves that W, V and T are invertible: for instance,
Wx=0 implies 60x=UWx=0. It avoids relying on a numerical rank calculation
or even on a particular symbolic determinant algorithm.

For a second form of the rank certificate, exact integer Bareiss elimination
with the stored row swaps and pivots gives

\[
\det W=-491520,\qquad
\det V=\frac{-491520}{60000}=-\frac{1024}{125},\qquad
\boxed{\det T=-\frac1{8000}}.
\tag{F14}
\]

All Bareiss divisions are checked to be exact. The searcher's separate
rational Gaussian elimination and an independent SymPy reconstruction give
the same determinants. The ordinary replay uses only the Python standard
library and checks (F13) in both multiplication orders.

## 6. Theorem and proof of the obstruction

**Theorem.** With \(\mathcal W,\mathcal C,\mathcal M\) defined in Section 2 and the prescribed
matching factors (F4)–(F7),

\[
\operatorname{span}_{\mathbb R}(\mathcal C\cup\mathcal M)
=\mathcal W,\qquad \dim\mathcal W=36.
\tag{F15}
\]

Consequently any real linear test subspace of W containing every member of
this specified family equals W. A test-space restriction that retains all
these functions cannot decrease the matrix order from 36.

**Proof.** The constant functions (F2) belong to W. For every admissible
matching parameter choice with delta positive, the direct multiplication
(F8) shows that the first complex top vector uses only the allowed
frequencies 0 and −2, the second uses only 0 and +2, and the bottom vector is
constant. For delta zero, (F9) gives the same containment directly. This
proves \(\operatorname{span}_{\mathbb R}(\mathcal C\cup\mathcal M)\subseteq\mathcal W\) over the full parameter
domains and both branches.

Conversely, the four functions (F11) belong to C. Each of the remaining 32
functions in the table is admissible for \(\mathcal M\), for the reasons given
after the table. Equations (F13) or (F14) prove that their 36 coefficient
columns in basis (F1) are linearly independent over the real numbers. Thus
their span has dimension 36. Since W has dimension 36, it equals W.
Finally, any linear subspace containing the entire specified family contains
these 36 independent functions and therefore cannot be a proper subspace
of W. This proves the theorem. □

This is a linear-span obstruction. It does **not** say that nonnegativity on
the 36 individual listed functions is equivalent to positivity of the
36 by 36 moment pencil: positivity on a collection of vectors need not imply
positivity on all their linear combinations. The conclusion is about which
linear test space can contain the original family.

## 7. Preservation of the published cone result

No restricted pencil is constructed because the first requested stopping
condition has been met. The existing lift remains the original homogeneous
36 by 36 pencil with 24 auxiliaries and 23 prescribed outputs. Its two
implications, the degree-20 root-based cone definition, the finite evaluation
witnesses, and its singular-boundary proof are unchanged.

In particular, the range conditions for all four B(r_j), the pseudoinverse
inequality, and the completed-square argument remain those of Section 8 of
the pinned proof. This span calculation does not replace or approximate
them. The z values used in the finite witness already satisfy the matching
hypothesis; no limit of matching factors, positive-definite-only
approximation, measure assumption, or closure operation enters this result.

The obstruction is specific to restricting a **real linear test space while
retaining every function from the fixed matching construction and the
constant family**. It makes no assertion about other matching factors,
smaller sufficient certificate families, other proof methods, different
lifts, or global matrix-size or auxiliary-count lower bounds. No novelty
claim is made.

## 8. Replay, independent checks, and the bounded stopping record

From the extracted handoff directory, run:

```bash
python verify_span.py
```

This requires Python 3.10 or later, with no third-party package, network
access, solver, or floating-point arithmetic. Do not run with Python's `-O`
option: the program rejects disabled assertions. The default run is read
only. To write a fresh report, use `--report new_report.json`.

The checker verifies the pinned source hashes and coordinate ordering,
reconstructs every selected test by **formal Laurent multiplication of the
source matrices**, checks exact parameter admissibility and branch values,
checks the rotation's unitarity and normalization and its two row identities,
checks the allowed frequency supports, and compares every saved matrix
entry. It then verifies the integer inverse identity (F13) and the exact
determinants (F14). It does not execute the search or any previous pencil
coefficient/kernel audit.

There are three independently implemented paths to the finite conclusion:

| Path | Construction | Exact rank check |
|---|---|---|
| Original bounded search | Gaussian-rational compact expansion | Rational Gaussian elimination |
| Primary replay | Separate Gaussian-rational field and formal Laurent matrix multiplication | Integer inverse identity and fraction-free Bareiss elimination |
| Independent review | SymPy expressions for the actual source R_0 over Q(sqrt(2),i), with no imports of local code | SymPy exact domain elimination |

The optional last replay is

```bash
python independent/check_fixed_witness_sympy.py
```

It requires SymPy; it was checked with SymPy 1.14.0. It only checks the fixed
witness. The primary verifier remains dependency free.

The search grid was declared and written before rank computation. It had
at most 1,881 matching candidates, in fixed angle/rho/u/D_sym order. The
search stopped immediately after the 638th candidate, which completed rank
36. Remaining candidates were not examined; the witness was not optimized
after discovery. `independent/declared_grid.json`, `search_log.json`, and the
retained search program document the bound and stopping rule. Replay of
the search is unnecessary for checking the result.

All reported ranks and identities are exact. The scripts check finite
algebraic claims and are not a proof-assistant formalization of every logical
step. The full-domain containment and implication from an invertible witness
to the span theorem are the elementary written arguments in Section 6.

The investigation is complete. No publication, repository mutation, external
communication, or scheduled follow-up was made.
