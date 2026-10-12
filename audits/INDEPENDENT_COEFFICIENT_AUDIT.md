# Independent coefficient reconstruction and bounded redundancy audit

Version: 1.0.0. This audit concerns the fixed 36-dimensional test space in
the task. It does not search other test spaces or claim a minimal lift.

The source read for the reconstruction is the supplied copy of
`build/sections/{introduction,angular,slices,matching,lift}.tex`, pinned to
OpenAI math commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`. The primary
generator was not imported or used to reconstruct the matrices. Only after
the independent reconstruction was complete was its exact output compared
with the primary `coefficients.json`.

## 1. Source and coordinate conventions

The 20-dimensional slice vector has five consecutive 4-dimensional
blocks, in the source order

\[
(\operatorname{Re}a,\operatorname{Im}a,
  \operatorname{Re}b,\operatorname{Im}b,u).
\]

Its slice pencil is

\[
 L_\theta=\begin{pmatrix}I_4\otimes\mathcal X&\mathcal B_\theta\\
 \mathcal B_\theta^T&\mathcal Z\end{pmatrix},\qquad
 \mathcal B_\theta=\begin{pmatrix}B(r_1)\\B(r_2)\\B(r_3)\\B(r_4)\end{pmatrix},
 \quad B(r)=\begin{pmatrix}0&r^T\\-r&0\end{pmatrix}.
\]

Here the columns \(r_1,\ldots,r_4\) are
\(\operatorname{Re}F_1,\operatorname{Im}F_1,
\operatorname{Re}F_2,\operatorname{Im}F_2\), with

\[
 U=\frac1{\sqrt2}\begin{pmatrix}-w^{-1}&-w\\
 s w^{-2}&s w^2\\c w^{-2}&-c w^2\end{pmatrix},\quad
 V=\frac1{\sqrt2}\begin{pmatrix}0&0\\i&i\\1&1\end{pmatrix},\quad
 F=\eta U+\xi V,
\]

where \(w=e^{i\theta}\), \(c=\cos\theta\), and \(s=\sin\theta\).
These formulas, including the signs in \(B\), are transcribed directly
in `source_pencil_parts()`.

The 36 test coefficients are ordered as

\[
 A_{\rm re}, A_{\rm im}, B_{\rm re}, B_{\rm im},
 C_{\rm re}, C_{\rm im}, D_{\rm re}, D_{\rm im}, u,
\]

with four coordinates per group, each in coordinate order 0,1,2,3.
They represent

\[
a=A+B w^{-2},\qquad b=C+D w^2.
\]

For an explicit real basis, let \(e_i\) be the four standard coordinate
vectors and put \(c_2=\cos2\theta\), \(s_2=\sin2\theta\). The nine groups
of basis vectors are the following five-block columns, for \(0\le i<4\):

| Group | Basis vector in the five source blocks |
|---|---|
| \(A_{\rm re},i\) | \((e_i,0,0,0,0)\) |
| \(A_{\rm im},i\) | \((0,e_i,0,0,0)\) |
| \(B_{\rm re},i\) | \((c_2e_i,-s_2e_i,0,0,0)\) |
| \(B_{\rm im},i\) | \((s_2e_i,c_2e_i,0,0,0)\) |
| \(C_{\rm re},i\) | \((0,0,e_i,0,0)\) |
| \(C_{\rm im},i\) | \((0,0,0,e_i,0)\) |
| \(D_{\rm re},i\) | \((0,0,c_2e_i,s_2e_i,0)\) |
| \(D_{\rm im},i\) | \((0,0,-s_2e_i,c_2e_i,0)\) |
| \(u,i\) | \((0,0,0,0,e_i)\) |

The basis is linearly independent: the complex frequencies 0 and -2
are independent in the first complex sector, and 0 and +2 are independent
in the second. The bottom constant vectors are independent of the top
sectors. Thus its real dimension is exactly 36.

For this basis, the matrix entries are exactly

\[
 H_{ij}(\Lambda)=\Lambda(q_i^T L_\theta q_j).
\]

There is no extra factor of two in an off-diagonal entry. In the
quadratic form \(t^THt\), the two symmetric entries provide the factor
of two. The independent verifier expands both directed contractions
\(q_i^TLq_j\) and \(q_j^TLq_i\) and checks their equality.

## 2. Exact frequency calculation

For the paired first two top blocks, the quadratic contribution is
\(\operatorname{Re}(a^*\mathcal X a)\); the second pair contributes
\(\operatorname{Re}(b^*\mathcal X b)\). Each uses only frequencies
\(0,\pm2\). Polarization yields the same support for every bilinear
entry. This pairing is why terms of frequency four cancel. The bottom
test is constant, so its \(\mathcal Z\) contribution is constant.

For an independent check of the mixed terms, the eight 4-by-4 top-to-bottom
blocks, before applying the functional, are

\[
\begin{split}
&B(\operatorname{Re}F_1),\quad B(\operatorname{Im}F_1),\quad
B(\operatorname{Re}(w^2F_1)),\quad B(\operatorname{Im}(w^2F_1)),\\
&B(\operatorname{Re}F_2),\quad B(\operatorname{Im}F_2),\quad
B(\operatorname{Re}(w^{-2}F_2)),\quad B(\operatorname{Im}(w^{-2}F_2)).
\end{split}
\]

Their \(\eta\) vectors have only harmonics 1 and 3, since the unshifted
columns of \(U\) have those harmonics and

\[
w^2 U_1=\frac1{\sqrt2}(-w,s,c)^T,\qquad
w^{-2} U_2=\frac1{\sqrt2}(-w^{-1},s,-c)^T.
\]

Their \(\xi\) vectors have only harmonics 0 and 2. The exact expansion
therefore lies in the 47-dimensional space with the following independent
coordinate functions:

| Slice variable | Real trigonometric basis | Dimension |
|---|---|---:|
| each of 10 upper-triangular \(\mathcal X_{ij}\) | \(1,\cos2\theta,\sin2\theta\) | 30 |
| each of 10 upper-triangular \(\mathcal Z_{ij}\) | \(1\) | 10 |
| \(\eta\) | \(\cos\theta,\sin\theta,\cos3\theta,\sin3\theta\) | 4 |
| \(\xi\) | \(1,\cos2\theta,\sin2\theta\) | 3 |

The verifier establishes both containment and that every displayed
frequency occurs. All calculations take place in
\(\mathbb Q(\sqrt2,i)[w,w^{-1}]\). Real trigonometric conversion is exact:
if the coefficient of \(w^k\), \(k>0\), is \(\alpha+i\beta\), the
cosine and sine coefficients are \(2\alpha\) and \(-2\beta\),
respectively. Conjugate symmetry is checked before this conversion.

## 3. Matching-certificate containment

The source matching construction has

\[
 M=\begin{pmatrix}a&b\\-\overline b&\overline a\end{pmatrix},
\]

where both \(a\) and \(\overline b\) have support \(\{0,-2\}\).
Its first column has that support; its second has support \(\{0,+2\}\).
This includes the source's exceptional \(\delta=0\) choice \(M=I_2\).
Consequently the first column of \(C=R_0M\) has the form
\(c_0+c_{-2}w^{-2}\) and its second has the form
\(d_0+d_2w^2\).

For every real \(u\) and every real symmetric \(D\), applying the fixed
real linear map \(-DJ_u^T\) to these columns preserves those supports
and commutes with real and imaginary parts. The source's entire family
of matching tests thus lies in the displayed 36-dimensional test space.
No positive-semidefiniteness condition on \(D\) is used here.

This containment is the needed mathematical sufficiency observation;
frequency counting alone would not establish the exact cone projection.
The closed-cone and singular-boundary argument is stated in the main
proof, and is outside the narrower scope of this coefficient program.

## 4. Outputs and the single bounded redundancy audit

All matrix and certificate indices below are zero based. Moments are
ordered by upper-triangular lexicographic \(X\) coordinate, with
\((1,\cos2,\sin2)\) for each, then the ten constant \(Z\) coordinates,
then \((\eta\cos1,\eta\sin1,\eta\cos3,\eta\sin3)\), and finally
\((\xi,\xi\cos2,\xi\sin2)\).

The output moment columns are

\[
(0,3,6,9,12,15,18,21,24,27,30,31,32,33,34,35,36,37,38,39,40,41,44).
\]

They are distinct coordinates, and their coordinate-selector minor is
\(I_{23}\). The remaining columns are the 24 auxiliary coordinates.

### Common vector kernel of the full and output coefficient families

Let \(H_k\) denote coefficient matrix \(k\). Exact addition gives

\[
H_0+H_{12}+H_{21}+H_{27}+H_{30}+H_{34}+H_{37}+H_{39}=I_{36}.
\]

This is the functional choice \(\Lambda(\mathcal X)=I_4\),
\(\Lambda(\mathcal Z)=I_4\), with every other moment zero. If a vector
lies in the common kernel of all coefficient matrices, the displayed
identity forces that vector to be zero. The exact common kernel is
therefore \(\{0\}\), with dimension zero and empty basis.
Every matrix in this displayed identity is an output coefficient, so
the output coefficient family by itself also has zero common vector
kernel.

### Auxiliary map with outputs fixed

The auxiliary coefficient map is

\[
\mathbb R^{24}\longrightarrow\operatorname{Sym}^{36},\qquad
v\longmapsto\sum_{k\in\mathrm{aux}}v_kH_k.
\]

Take the following 24 matrix-entry rows, in the corresponding auxiliary
order. For each \(0\le i\le j<4\), take

\[
(A_{{\rm re},i},B_{{\rm re},j}),\qquad
(A_{{\rm re},i},B_{{\rm im},j}).
\]

These select the \(X_{ij}\cos2\) and \(X_{ij}\sin2\) coefficients.
Append the four rows

\[
(A_{{\rm im},0},u_2),\quad (A_{{\rm re},0},u_2),\quad
(B_{{\rm im},0},u_2),\quad (B_{{\rm im},0},u_3).
\]

After restricting columns to auxiliaries, the resulting minor is exactly

\[
\operatorname{diag}\left(
\underbrace{1,\ldots,1}_{20},
\frac{\sqrt2}{4},\frac{\sqrt2}{4},
\frac{\sqrt2}{2},\frac{\sqrt2}{2}\right),
\]

with determinant \(1/16\). Hence the auxiliary map has exact rank 24
and zero kernel. No auxiliary direction is redundant in this sense.

### Common vector kernel of the auxiliary coefficient family

The parameter-map kernel above consists of vectors in \(\mathbb R^{24}\).
It must not be confused with the intersection of the kernels of the
auxiliary coefficient matrices, whose vectors lie in \(\mathbb R^{36}\).
For the auxiliary coefficient matrices, this common vector kernel is
exactly

\[
 \bigcap_{k\in\mathrm{aux}}\ker H_k=\operatorname{span}\{e_{33}\}.
\]

Here \(e_{33}\) is the zero-based coordinate for \(u_1\) in the 36-test
ordering. To certify the statement within the same bounded audit, set

\[
 P=H_1+H_{13}+H_{22}+H_{28}.
\]

These are the four diagonal \(X\cos2\) coefficients. Direct exact
multiplication gives \(P^2=\operatorname{diag}(I_{32},0_4)\), so any
common-kernel vector has zero first 32 coordinates. For \(H_{45}\),
the \(\xi\cos2\) coefficient, the submatrix in rows \((11,12,8)\)
and columns \((32,34,35)\) is

\[
 \operatorname{diag}\left(-\frac{\sqrt2}{2},
 \frac{\sqrt2}{2},\frac{\sqrt2}{2}\right),
\]

with determinant \(-\sqrt2/4\). Thus only coordinate 33 can survive.
The program checks that every auxiliary coefficient annihilates
\(e_{33}\), proving equality with the stated one-dimensional space.

This does not permit deletion of coordinate 33 from the full pencil:
the output coefficient \(H_{34}\), the constant \(Z_{22}\) moment in
one-based matrix-entry notation, maps \(e_{33}\) to itself. The full
coefficient family has zero common kernel, as certified above.

### Full coefficient map

For completeness, a 47-by-47 minor of the full map is also specified
and verified. For every upper-triangular \((i,j)\), its first 30 rows
are \((i,j),(i,8+j),(i,12+j)\), in that order for each pair. The next
ten rows are \((32+i,32+j)\). The last seven rows are

\[
(8,33),(12,33),(4,34),(0,34),(20,34),(12,34),(12,35).
\]

Its columns are all 47 moments in canonical order. Its exact determinant
is \(\sqrt2/32\). Thus the full coefficient map is injective, of rank
47. In particular the bilinear test products span the entire claimed
47-dimensional function space, rather than a proper subspace.

This is one bounded audit of the stated pencil. It removes nothing,
because the full coefficient common vector kernel and the auxiliary
parameter-map kernel are zero. The auxiliary coefficient common vector
kernel described above is not a removable full-pencil direction.
It establishes no lower bound on the size of arbitrary alternative
lifts and excludes no other type of semidefinite reformulation.

## 5. Replay and results

Run:

```sh
python3 audits/independent_coefficients.py --primary-coefficients coefficients.json
```

The program writes:

- `independent_coefficients.json`: all 47 exact sparse upper-triangular
  coefficient matrices, using rational pairs for \(a+b\sqrt2\);
- `redundancy_certificate.json`: exact minor rows, columns, determinants,
  common-kernel identity, and rank/kernel statements;
- `independent_verification_log.json`: the independent expansion result;
- `independent_comparison.json`: exact comparison with the primary matrix
  coefficient file.

The completed run passed the following checks:

- all 28,512 directed contractions \(q_i^T L^{(v)}q_j\), corresponding
  to 36 squared ordered test pairs and 22 slice-coordinate matrices;
- exact reality, symmetric polarization, and frequency supports,
  including complete lower-triangle coefficient-dictionary comparisons
  that reject every unexpected mode;
- all 47 matrices and all 60,912 full matrix scalar entries against the
  independently implemented compact primary generator;
- the 298 nonzero upper-triangular scalar entries;
- output rank 23, full coefficient-map rank 47, auxiliary parameter-map
  rank 24, full/output common vector kernel zero, and auxiliary
  coefficient common vector kernel \(\operatorname{span}\{e_{33}\}\),
  with the certificates above.

These are symbolic algebra checks. They do not use numerical tolerances,
angular samples, random tests, or a representing-measure hypothesis.
They do not by themselves establish the lift theorem; the associated
mathematical proof supplies that conclusion. Originality is not assessed
by this audit.
