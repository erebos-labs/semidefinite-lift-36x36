# Independent audit of the matching test family

## Scope

This note uses the authoritative `source/THEOREM_PROOF.md` and
`source/COORDINATES.md` from the stated erebos pin
`2d2c170206ab937bcf56ff3290066cc3b2dfdd99`. It concerns only the real linear
span of the tests constructed in equation (26). It does not repeat the old
coefficient redundancy audit, search for other lifts, or change the completed
36/24 construction.

The source fixes a real 36-coordinate order

\[
(\Re A,\Im A,\Re B,\Im B,\Re C,\Im C,\Re D_f,\Im D_f,u),
\]

with four coordinates per group and
`a = A + B exp(-2 i theta)`, `b = C + D_f exp(+2 i theta)`.
Here `D_f` is a complex Fourier coefficient. The real symmetric matrix `D`
used to construct a test is a different object.

## 1. Precisely admissible inputs

A constructed matching test is determined by

\[
u\in\mathbb R^4,\qquad D=D^T\in\mathbb R^{4\times4},\qquad
z=\tau(c,s,0)+\rho(0,0,1),
\]

where `tau > 0`, `c^2+s^2=1`, and `rho` is real. Thus `(z1,z2)` is nonzero.
There is no definiteness condition on `D` or on `S=J_u D J_u^T`.
The source's matrix `S` is determined by `u,D`; it is not an additional
independently selectable input. Likewise `-D J_u^T R0` cannot be replaced
by a free complex matrix when constructing a span witness.

Set `w0=c+i s`, define `R0=tau U(alpha)+rho V0`, and use exactly the source's
`g,h,delta` derived from `N=S R0`. The matching factor is the explicit
piecewise factor in Section 6 of the source. Equation (26) then has complex
top pairs

\[
(a(\theta),b(\theta))=-D J_u^T R_0M(\theta)
\]

and bottom block `u`.

The `z`-axis is not part of this explicit parameter family. The old cone
proof extends a scalar inequality to that axis; it does not supply an
axis-valued matching factor by continuity. A finite span certificate must
use admissible nonaxis parameters or separately prove that a proposed extra
test is in the span of admissible tests.

## 2. Exact 36-coordinate formula

Put

\[
K=-D J_u^T,\qquad (p,q)=K R_0,\qquad t=w_0^2.
\]

For the `delta > 0` branch define

\[
r=\frac{|g|^2}{\delta},\qquad
s_0=\frac{|h|^2}{\delta},\qquad v=\frac{gh}{\delta}.
\]

The four complex coefficient columns are

\[
\begin{aligned}
A&=rp-vq,\\
B&=t(s_0p+vq),\\
C&=\overline v\,p+rq,\\
D_f&=\overline t(-\overline v\,p+s_0q).
\end{aligned}
\tag{A}
\]

Indeed, the first column of `M` is `(a0,-conjugate(b0))^T`, where
`a0=r+s0 t w^-2` and `conjugate(b0)=v(1-t w^-2)`. Its second column is
their required conjugate-frequency counterpart. Multiplying gives (A)
without sampling angles or inferring Fourier coefficients numerically.

For the **separate** `delta = 0` branch, the source prescribes `M=I`, so

\[
A=p,\qquad B=0,\qquad C=q,\qquad D_f=0.
\tag{B}
\]

There is no division or limiting interpretation in (B). In either branch,
take real and imaginary parts in the exact nine-group coordinate order
displayed above. In particular, the imaginary part of the coefficient of
`w^-2` is used directly: its sign is already encoded by the source's real
basis `(sin(2theta)e_i,cos(2theta)e_i,0,0,0)`.

### Rational normalization used in a finite witness

If `c,s,tau,rho,u,D` are rational, with `c^2+s^2=1`, define
`Rhat=sqrt(2) R0`. This matrix is Gaussian rational. Derive
`ghat,hhat` from `Nhat=S Rhat` using the same formulas. Then

\[
\widehat g=\sqrt2g,\quad\widehat h=\sqrt2h,\quad
\widehat\delta=2\delta.
\]

The three ratios `r,s0,v` are unchanged. Using
`(phat,qhat)=K Rhat` in (A) or (B) produces `sqrt(2)` times each of the
32 true top coordinates. Keeping the bottom `u` unchanged produces a
fully rational vector obtained from the true coefficient vector by

\[
T=\operatorname{diag}(\sqrt2 I_{32},I_4).
\]

This is an invertible transformation. It preserves rank, but its determinant
is `2^16`, which must be included when translating a 36 by 36 determinant
back to the source coordinates. The test parameters and their true source
coordinates should not be confused with these rationally scaled coordinates.

## 3. Useful exact inputs for both branches

The stated bounded grid `w0 in {1,i,(3+4i)/5}`, `tau=1`, and
`rho in {0,1,-1}` has valid nonaxis `z` throughout. Basis vectors and sums
for `u`, and arbitrary explicitly symmetric rational matrices for `D`,
are allowed. Exact computation must decide which denominator branch each
individual input occupies.

For example, `D=0`, `u=e_i`, `z=(1,0,0)` yields `S=0`, `delta=0`, and
the pure bottom test `(0,0,0,0,e_i)`. A nonzero indefinite zero-denominator
example is `u=e_0`, `D=diag(2,1,-1,0)`, `z=(1,0,0)`.

Changing the real symmetric matrix from `D` to `lambda D`, for real
`lambda != 0`, leaves the rotation ratios unchanged: `g,h` both scale by
`lambda`, while `delta` scales by `lambda^2`. Its top 32 coordinates scale
by `lambda`, and its bottom coordinates remain `u`. This observation can
be used algebraically, but a witness should still record the actual
admissible parameters for every test it includes.

## 4. The constant family and the exact scope of the conclusion

The four selected constant tests are independently admissible members of
`C={(v,0,0,0,0):v in R^4}`, the constant family specified in the request.
They are included explicitly. No assertion that they must themselves be
single equation-(26) tests is needed.

Every source matching test is in the fixed 36-dimensional space `W`, as
shown by formulas (A), (B). A finite rank-36 matrix containing four constants
and 32 admissible matching tests proves that the real span of the augmented
family `C union M` is `W`. Thus no proper real linear subspace of `W` can
retain all members of this specified family.

It would not prove that the 36 by 36 lift is globally minimal, that 24
auxiliaries are minimal, or that every individual matching test is logically
necessary for the cone argument. A different choice of certificates or a
different proof could use different tests. This result concerns a fixed
vector family, not all lifts or all possible proof reorganizations.

Also, spanning by finitely many test vectors does not make positivity on
those finitely many vectors equivalent to positive semidefiniteness of the
whole moment matrix. For example, two coordinate vectors span `R^2`, but
the symmetric matrix `[[1,-2],[-2,1]]` is positive on each coordinate
vector and negative on their sum. No such positivity equivalence is
claimed from the rank computation.

## 5. Independent replay of the supplied fixed witness

The subsequently supplied `independent/span_witness.json` was independently
reconstructed by `check_fixed_witness_sympy.py`, using the actual source
coefficients in `Q(sqrt(2),i)` and formulas (A), (B). The checker imports no
local generator, search implementation, or other verifier and performs no
search. All 1,296 entries agreed with the saved normalized columns; the
row-major and column-major copies also agreed.

The composition is four constant first-block tests and 32 constructed
matching tests. Of the matching tests, eight use `delta=0` and 24 use
`delta>0`. All have exactly symmetric `D`, `tau>0`, unit-circle angle
parameters, and nonaxis `z`; the recorded `g_hat,h_hat,delta,delta_hat`
were checked against independent source calculations.

SymPy exact rational elimination gave normalized determinant `-1024/125`.
The row-scaling determinant is `2^16=65536`, so the determinant in the
actual source coordinate basis is **`-1/8000`**, which is nonzero. The
supplied augmented proof-test family therefore has real span dimension
**36**. This is the bounded follow-up's requested result.

The independently reviewed witness SHA-256 is
`481e5c5ef4ab866624ac77ca737bd0fb7dab73af7a300dc417f3be27f068d327`.
The replay record is `independent/sympy_witness_review.json`.

Replay with SymPy 1.14.0 available:

```sh
python independent/check_fixed_witness_sympy.py
```

No additional candidates, alternative lift constructions, or coefficient
redundancy analyses were investigated during this fixed-witness review.
