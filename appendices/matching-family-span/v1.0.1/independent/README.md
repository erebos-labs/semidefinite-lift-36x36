# Finite exact witness for the matching-test span

## Result

The search found **36 independent admissible tests** in the specified
36-dimensional real test space: the four constant tests
\(e_i/\sqrt2\), \(0\le i<4\), supported in the first four-vector block,
and **32 tests from the authoritative matching construction**.

Of the matching tests selected, **8 use the exact \(\delta=0\) branch**
\(M=I_2\), and **24 use the \(\delta>0\) formula**. No limit between
branches is used.

The ordered coefficient matrix has exact determinant

\[
\det T=-\frac1{8000}\ne0.
\]

After multiplying its first 32 rows by \(\sqrt2\), the matrix is rational
and its exact determinant is

\[
\det\bigl(\operatorname{diag}(\sqrt2I_{32},I_4)T\bigr)
=-\frac{1024}{125}.
\]

The determinant scale factor is \((\sqrt2)^{32}=2^{16}\). Both matrices
refer to the 36-coordinate order in `source/COORDINATES.md`.

Every exhibited test belongs to the prescribed space by the support of
the authoritative matching formula. The nonzero determinant proves that
these admissible matching tests together with the specified constants
span the entire space. This statement concerns the span of these tests;
it gives no minimality conclusion for arbitrary semidefinite lifts.

## Authoritative provenance

The construction is transcribed from
`erebos-labs/semidefinite-lift-36x36`, commit
`2d2c170206ab937bcf56ff3290066cc3b2dfdd99`:

- `THEOREM_PROOF.md`, Section 3: ordered test basis and equation (6);
- Section 6: both branches of the matching construction, especially (22);
- Section 8: the admissible test (26), with \(S=J_uDJ_u^T\).

The old lift coefficient and redundancy checks are not rerun here.
`search_matching_span.py` uses a compact coefficient expansion of these
formulas. The root verifier independently reconstructs the selected tests
by source-level Laurent matrix multiplication; its record is separate.

## The declared bound and stopping point

The grid was fixed and communicated before execution, and is preserved in
`declared_grid.json`:

- \(\tau=1\);
- \(w_0\in\{1,i,(3+4i)/5\}\), in that order;
- \(\rho\in\{0,1,-1\}\), in that order;
- \(u\): the four standard basis vectors, the six pairwise sums in
  lexicographic order, and \((1,1,1,1)\);
- \(D\): zero, identity, four diagonal units, six symmetric off-diagonal
  units, six identity-plus-off-diagonal units, and
  \(\operatorname{diag}(1,-1,1,-1)\), in that order.

The nested order is \(w_0,\rho,u,D\), giving at most **1,881 matching
candidates**. Search stopped immediately at rank 36, after **638
candidates**, with last zero-based candidate index **637**. No remaining
candidates were examined. The search is complete and stopped; no grid
expansion, witness optimization, or alternative-lift search was performed.

Only angles \(w_0=1\) and \(w_0=i\), and values \(\rho=0,1\), occur
in the selected witness. A total of 280 zero-branch and 358 positive-branch
candidates were examined before the stopping point.

## Selected exact parameters

The first four columns are the first-block constants \(e_i/\sqrt2\),
\(i=0,1,2,3\). The table gives the remaining columns. All row, column,
and vector indices in this table are zero based. Let
\(E_{ii}=e_ie_i^T\), \(F_{ij}=e_ie_j^T+e_je_i^T\), and let
\(\mathbf1=(1,1,1,1)^T\). Every row has \(\tau=1\), and
\(\widehat\delta=2\delta\). Thus a zero in the final column specifies
the separate \(M=I_2\) branch.

| Column | \(w_0\) | \(\rho\) | \(u\) | Symmetric \(D\) | \(\widehat\delta\) |
|---:|---|---:|---|---|---:|
| 4 | 1 | 0 | \(e_0\) | 0 | 0 |
| 5 | 1 | 0 | \(e_0\) | \(I\) | 4 |
| 6 | 1 | 0 | \(e_0\) | \(E_{11}\) | 0 |
| 7 | 1 | 0 | \(e_0\) | \(E_{33}\) | 4 |
| 8 | 1 | 0 | \(e_0\) | \(F_{01}\) | 0 |
| 9 | 1 | 0 | \(e_0\) | \(F_{12}\) | 4 |
| 10 | 1 | 0 | \(e_0\) | \(F_{13}\) | 4 |
| 11 | 1 | 0 | \(e_0\) | \(F_{23}\) | 4 |
| 12 | 1 | 0 | \(e_0\) | \(I+F_{01}\) | 4 |
| 13 | 1 | 0 | \(e_0\) | \(I+F_{03}\) | 4 |
| 14 | 1 | 0 | \(e_0\) | \(I+F_{13}\) | 8 |
| 15 | 1 | 0 | \(e_0\) | \(I+F_{23}\) | 8 |
| 16 | 1 | 0 | \(e_1\) | 0 | 0 |
| 17 | 1 | 0 | \(e_2\) | 0 | 0 |
| 18 | 1 | 0 | \(e_3\) | 0 | 0 |
| 19 | 1 | 0 | \(e_3\) | \(I+F_{02}\) | 4 |
| 20 | 1 | 0 | \(e_0+e_1\) | \(I+F_{23}\) | 8 |
| 21 | 1 | 0 | \(e_0+e_2\) | \(I+F_{03}\) | 8 |
| 22 | 1 | 0 | \(e_0+e_3\) | \(I+F_{01}\) | 20 |
| 23 | 1 | 0 | \(e_0+e_3\) | \(I+F_{02}\) | 20 |
| 24 | 1 | 0 | \(e_2+e_3\) | \(I+F_{01}\) | 8 |
| 25 | 1 | 0 | \(e_2+e_3\) | \(I+F_{03}\) | 8 |
| 26 | 1 | 0 | \(\mathbf1\) | \(F_{03}\) | 24 |
| 27 | 1 | 0 | \(\mathbf1\) | \(I+F_{12}\) | 40 |
| 28 | 1 | 1 | \(e_0\) | \(I\) | 20 |
| 29 | 1 | 1 | \(e_0\) | \(F_{02}\) | 0 |
| 30 | 1 | 1 | \(e_0\) | \(F_{12}\) | 4 |
| 31 | 1 | 1 | \(e_0\) | \(F_{23}\) | 4 |
| 32 | \(i\) | 0 | \(e_0\) | \(I\) | 4 |
| 33 | \(i\) | 0 | \(e_0\) | \(F_{01}\) | 0 |
| 34 | \(i\) | 0 | \(e_0\) | \(F_{12}\) | 4 |
| 35 | \(i\) | 0 | \(e_0\) | \(F_{13}\) | 4 |

## Files and arithmetic

- `span_witness.json` contains every exact parameter, branch,
  \(\widehat\delta\), the rational scaled 36-by-36 matrix in row and
  column formats, coordinate labels, and both determinants.
- `search_log.json` records the completed stopping point and branch counts.
- `declared_grid.json` specifies the fixed candidate grid.
- `search_matching_span.py` preserves the stopped construction procedure.

The construction uses Gaussian rationals for \(\widehat R_0=\sqrt2R_0\).
Since \(\widehat g=\sqrt2g\) and \(\widehat h=\sqrt2h\), the normalized
rotation formula is unchanged, and
\(\widehat\delta=|\widehat g|^2+|\widehat h|^2=2\delta\).
The first 32 actual test coefficients equal their stored rational values
times \(\sqrt2/2\); the bottom four coefficients equal their stored
rational values. All arithmetic uses Python's standard-library
`Fraction`. No numerical tolerance, genericity assumption, or sampled
angular values enter the rank calculation.

The primary delivery should use the separate witness verifier to check
the saved 36 columns directly. Repeating or extending the search is not
needed for verification.
