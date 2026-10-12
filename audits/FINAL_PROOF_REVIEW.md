# Final independent review of the main proof

**Artifact reviewed:** `THEOREM_PROOF.md`, version 1.0.0.

**Outcome:** PASS for the mathematical construction and proof. No unresolved
mathematical gap or unsupported dimensional assertion was found. One cosmetic
transpose typo in equation (21) was reported for correction: the exponent
`,T` should be `T`. It does not change the argument or generated matrices.

This final review independently read the complete main proof and checked its
new compact formulas against the pinned source, its explicit real test basis,
the exact generator, and the recorded coefficient/rank certificates. It also
read the separate denominator-free determinant checker. It did not perform
an additional search for alternative lifts or repeat already completed
verification runs.

## 1. The theorem states the correct result

The displayed `Q`, `Phi`, degree-20 `p`, direction `e`, and root convention
match the seven requested pinned source files. The theorem has 23 output
coordinates, a homogeneous 36 by 36 real symmetric pencil, and 24 independent
real auxiliary coordinates. Its coefficients lie in the stated field
`Q(sqrt(2))`. Both singular `X` and zero polynomial roots remain included.

The note explicitly distinguishes a real mathematical proof supplemented by
exact symbolic verification from a formal proof-assistant theorem. It also
confines the no-redundancy result to this particular pencil and makes no
global minimality assertion.

## 2. Compact top blocks and all eight mixed columns are correct

The ordered basis has first complex pair

`a = A + B exp(-2 i theta)`

and second complex pair

`b = C + D exp(+2 i theta)`.

For the first pair, the cross matrix between the constant and oscillating
real coefficient groups is

\[
\begin{pmatrix}X_c&X_s\\-X_s&X_c\end{pmatrix}.
\]

For the second pair it is

\[
\begin{pmatrix}X_c&-X_s\\X_s&X_c\end{pmatrix}.
\]

The oscillating self-pairings give `X` since the squared sine and cosine terms
sum to one; their real/imaginary cross pairing is zero. Thus the displayed
`T_+` and `T_-` are correct, including every sign.

I manually re-expanded all eight source-ordered mixed columns. They are
exactly the real and imaginary parts of

\[
f_1,\quad w^2f_1,\quad f_2,\quad w^{-2}f_2,
\]

in that order, with each real part followed by its imaginary part. The useful
shifted expressions are

\[
\sqrt2 w^2f_1=\eta(-w,s,c)^T+\xi(0,iw^2,w^2)^T,
\]

\[
\sqrt2 w^{-2}f_2=\eta(-w^{-1},s,-c)^T+
\xi(0,iw^{-2},w^{-2})^T.
\]

These give rows 3, 4, 7, 8 of the main table directly. For the remaining
rows, the identities

\[
s\cos2\theta=(\sin3\theta-\sin\theta)/2,
\quad s\sin2\theta=(\cos\theta-\cos3\theta)/2,
\]

\[
c\cos2\theta=(\cos3\theta+\cos\theta)/2,
\quad c\sin2\theta=(\sin3\theta+\sin\theta)/2
\]

give precisely the displayed `h_c`, `h_s`, and output terms. Every table
entry includes the correct common factor `1/sqrt(2)`.

The skew matrix `B(r)` has the correct source sign. Its use in upper mixed
blocks and its transpose in lower mixed blocks makes the full pencil
symmetric. The polarization formula is correct, with no extra factor of two
in an off-diagonal stored coefficient.

## 3. The denominator-free determinant certificate is valid

For each fixed real angle, work in the polynomial ring over the real
coefficients in the independent slice variables. In the main proof's notation,

\[
A\widehat J=dI_{16},\qquad
\mathcal B^T\widehat J\mathcal B
=\Phi_{Y_\theta}(\operatorname{adj}\mathcal X).
\]

Multiplication by the displayed right matrix therefore gives

\[
L_\theta R=\begin{pmatrix}
dA&0\\d\mathcal B^T&d\mathcal Z-
\Phi_{Y_\theta}(\operatorname{adj}\mathcal X)
\end{pmatrix}.
\]

The two powers in the determinant argument agree:

\[
\det R=d^{20},\qquad
\det(dA)=d^{16}\det A=d^{16}d^4=d^{20}.
\]

Thus `d^20 det(L_theta) = d^20 p` holds as a polynomial identity. The
polynomial `d` is nonzero (its value at the identity is one), and this
polynomial ring is an integral domain. Cancellation is therefore legitimate
before evaluation. The result can then be evaluated at singular `X`; it
does not divide by a singular evaluated determinant and does not remove any
factor from the defining polynomial.

The determinant checker verifies the full block product for arbitrary four
real skew columns and checks the adjugate identities. Its dependence on the
separately verified angular Gram identity is stated accurately. The root
identity that follows is the degree-20 characteristic polynomial of the
real symmetric slice pencil, so it preserves all roots and multiplicities.

## 4. Matching and both lift implications are complete

The main proof's matching construction reproduces the independently audited
source identities and retains the stronger one-sided column supports. For
arbitrary real symmetric `D`, the constant real map `-D J_u^T` preserves
those supports and the realification order, so every required matching
certificate belongs to the 36-dimensional test space.

The positive-denominator and zero-denominator branches are handled separately.
The derivative argument uses only symmetry, not positivity, of `S`. Its use
of the first harmonic's value and derivative is valid. A matching factor may
depend on `u`, `D`, and `z`; positivity against every test vector permits
these choices independently.

The forward implication supplies explicit finite evaluation moments for each
cone point. The converse starts from an arbitrary functional, uses constant
top tests to prove `X` positive semidefinite, and derives the scalar
certificate with arbitrary symmetric `D`. The `z`-axis extension concerns
only a scalar polynomial inequality; it neither assumes continuity of
matching factors nor takes a closure of a projection.

For singular `X`, the proof correctly obtains all four separate `B_j` range
conditions using arbitrary real multiples of `v v^T` for kernel vectors.
The pseudoinverse choice then supplies the required matrix inequality. The
completed-square identity has the correct cross-term orientation and uses
the already established relations `X Xdagger B_j = B_j`. This establishes
actual slice positivity and hence the original closed-cone membership.

## 5. The bounded redundancy conclusions have exact support

The eight diagonal output coefficient matrices listed in equation (31) sum
to `I_36`. Their indices are correct in the documented moment ordering.
This proves that the common matrix kernel of the full coefficient family is
zero.

For the auxiliary map, the selected top-block entries isolate the twenty
`X_c`, `X_s` coordinates. The four mixed entries isolate `h_c`, `h_s`,
`k_c`, `k_s` with coefficients

\[
\sqrt2/4,\quad\sqrt2/4,\quad\sqrt2/2,\quad\sqrt2/2.
\]

The resulting diagonal minor has determinant `1/16`, so its claimed rank
24 and zero parameter kernel follow. The saved full-map certificate records
rank 47 and determinant `sqrt(2)/32` for its specified minor, consistent
with the main proof.

The auxiliary matrices' common **vector** kernel is correctly distinguished
from their **parameter** kernel. Their top part is killed by requiring the
`X_c = I` combination to annihilate it. The remaining `k_c` blocks force
`u_0 = u_2 = u_3 = 0`, leaving only the zero-based coordinate 33. Every
auxiliary mixed column has zero first spatial component, so all auxiliary
matrices annihilate that last direction. Output `Z` acts on it; it cannot
be removed as a common kernel of the complete pencil.

## 6. Verification claims and review limits

The saved primary verification result reports all seven pinned source hashes,
47 size-36 matrices, 298 nonzero upper-triangle coefficients, and completed
independent comparison of 60,912 full matrix entries. These counts agree
with the dimensions and the main proof. The proof correctly treats the
exact singular examples as regressions rather than as universal evidence.

The completed exact matching checker supplies 42 assertions and is separate
from the primary generator. The written proof supplies the logical steps
outside those finite identities. This division of responsibility is stated
without presenting numerical tests as a proof.

Historical novelty and current-repository search claims are covered by the
separate prior-work audit; this review did not repeat that search. The main
proof's express qualifications about novelty, formalization, numerical
conditioning, and global minimality are appropriate.

## Resolution of the notation correction

The primary author corrected equation (21) to `\widetilde{\mathscr C}^{T}` before packaging v1.0.0. No coefficient, proof implication, or verification result changed.
