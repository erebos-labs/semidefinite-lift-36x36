# Coordinate ordering and coefficient format — v1.0.0

All machine indices below are zero based. Output X/Z entries use mathematical matrix indices 1 through 4; test vector coordinates use 0 through 3 as in the source u=(u0,u′).

## Test coordinates

The nine groups, each with four coordinates 0,1,2,3, are A.real, A.imag, B.real, B.imag, C.real, C.imag, D.real, D.imag, u. Here a=A+B exp(-2i theta), b=C+D exp(+2i theta). Realification is (Re a, Im a, Re b, Im b, u). The complete real basis is given in Section 3 of THEOREM_PROOF.md.

| Test index | Coefficient |
|---|---|
| 0 | `A.real.0` |
| 1 | `A.real.1` |
| 2 | `A.real.2` |
| 3 | `A.real.3` |
| 4 | `A.imag.0` |
| 5 | `A.imag.1` |
| 6 | `A.imag.2` |
| 7 | `A.imag.3` |
| 8 | `B.real.0` |
| 9 | `B.real.1` |
| 10 | `B.real.2` |
| 11 | `B.real.3` |
| 12 | `B.imag.0` |
| 13 | `B.imag.1` |
| 14 | `B.imag.2` |
| 15 | `B.imag.3` |
| 16 | `C.real.0` |
| 17 | `C.real.1` |
| 18 | `C.real.2` |
| 19 | `C.real.3` |
| 20 | `C.imag.0` |
| 21 | `C.imag.1` |
| 22 | `C.imag.2` |
| 23 | `C.imag.3` |
| 24 | `D.real.0` |
| 25 | `D.real.1` |
| 26 | `D.real.2` |
| 27 | `D.real.3` |
| 28 | `D.imag.0` |
| 29 | `D.imag.1` |
| 30 | `D.imag.2` |
| 31 | `D.imag.3` |
| 32 | `u.0` |
| 33 | `u.1` |
| 34 | `u.2` |
| 35 | `u.3` |

## Moment, output and auxiliary coordinates

Each label is a real functional value, for example X12.cos2 is Lambda(mathcal X12 cos(2 theta)). It is not a complex Fourier coefficient.

| Moment index | Label | Role and zero-based index |
|---|---|---|
| 0 | `X11.1` | output 0 |
| 1 | `X11.cos2` | auxiliary 0 |
| 2 | `X11.sin2` | auxiliary 1 |
| 3 | `X12.1` | output 1 |
| 4 | `X12.cos2` | auxiliary 2 |
| 5 | `X12.sin2` | auxiliary 3 |
| 6 | `X13.1` | output 2 |
| 7 | `X13.cos2` | auxiliary 4 |
| 8 | `X13.sin2` | auxiliary 5 |
| 9 | `X14.1` | output 3 |
| 10 | `X14.cos2` | auxiliary 6 |
| 11 | `X14.sin2` | auxiliary 7 |
| 12 | `X22.1` | output 4 |
| 13 | `X22.cos2` | auxiliary 8 |
| 14 | `X22.sin2` | auxiliary 9 |
| 15 | `X23.1` | output 5 |
| 16 | `X23.cos2` | auxiliary 10 |
| 17 | `X23.sin2` | auxiliary 11 |
| 18 | `X24.1` | output 6 |
| 19 | `X24.cos2` | auxiliary 12 |
| 20 | `X24.sin2` | auxiliary 13 |
| 21 | `X33.1` | output 7 |
| 22 | `X33.cos2` | auxiliary 14 |
| 23 | `X33.sin2` | auxiliary 15 |
| 24 | `X34.1` | output 8 |
| 25 | `X34.cos2` | auxiliary 16 |
| 26 | `X34.sin2` | auxiliary 17 |
| 27 | `X44.1` | output 9 |
| 28 | `X44.cos2` | auxiliary 18 |
| 29 | `X44.sin2` | auxiliary 19 |
| 30 | `Z11.1` | output 10 |
| 31 | `Z12.1` | output 11 |
| 32 | `Z13.1` | output 12 |
| 33 | `Z14.1` | output 13 |
| 34 | `Z22.1` | output 14 |
| 35 | `Z23.1` | output 15 |
| 36 | `Z24.1` | output 16 |
| 37 | `Z33.1` | output 17 |
| 38 | `Z34.1` | output 18 |
| 39 | `Z44.1` | output 19 |
| 40 | `eta.cos1` | output 20 |
| 41 | `eta.sin1` | output 21 |
| 42 | `eta.cos3` | auxiliary 20 |
| 43 | `eta.sin3` | auxiliary 21 |
| 44 | `xi.1` | output 22 |
| 45 | `xi.cos2` | auxiliary 22 |
| 46 | `xi.sin2` | auxiliary 23 |

## Exact scalar and matrix encoding

`coefficients.json` stores one array for each of the 47 moments in the table. An entry `[row,col,[a_num,a_den,b_num,b_den]]` means matrix value a+b sqrt(2), where a=a_num/a_den and b=b_num/b_den. Denominators are positive and fractions are reduced.

Only nonzero upper-triangular entries are stored. Every omitted entry is zero. Reflect off-diagonal entries into the lower triangle unchanged: no additional factor two. The coefficient of an off-diagonal output matrix entry Xij or Zij already includes both symmetric placements in that 4 by 4 matrix.

The pencil is the linear sum of the 47 stored matrices multiplied by the corresponding moment values. Select the 23 output matrices by `output_indices` and the 24 auxiliary matrices by `auxiliary_indices`. There is no affine constant term.

The primary matrices use scalar values 0, ±1, ±sqrt(2)/2 and ±sqrt(2)/4. General exact arithmetic in generator.py uses Q(sqrt(2)) and intentionally rejects floats.

## Independent file format

`audits/independent_coefficients.json` uses a deliberately separate encoding: each sparse entry is `[row,col,a,b]` with rational strings a,b for a+b sqrt(2). The independent reconstruction converts and compares every coefficient against the primary file. The different file layouts are documented and are not differences in the matrices.

## Rank certificate indices

`certificates/rank_kernel.json` uses the primary moment ordering and zero-based matrix-entry rows. `audits/redundancy_certificate.json` records independently chosen row orders for its minors; use the row order stored with each determinant. Both full-map minors have determinant sqrt(2)/32, and both auxiliary minors have determinant 1/16.

The common vector kernel of the full matrix coefficient family, the common vector kernel of the auxiliary subfamily, and the kernel of the map from auxiliary parameters to matrices are distinct objects; their dimensions are respectively 0, 1, and 0.
