# Independent model audit — public abridgment

This is a sanitized, abridged copy of the completed private mathematical audit, not a raw original log or a human peer review. Acquisition/account identifiers, machine paths and operational notes are omitted. Mathematical reasoning and reported results in the retained sections are unchanged. No specific theorem gap was found. The guard-only verifier correction was subsequently tested; public release controls rerun it. Nonspectrahedrality is inherited from the companion and was not independently audited.

## First-principles mathematical review

### Restricted test space, moments, and pencil

The first complex pair is A+B w^-2 and the second is C+D w^2. Each has two independent complex four-vector coefficients, hence 16 real coordinates. Adding constant real u gives 36. Realification is exactly (Re a, Im a, Re b, Im b, u), with the stated signs of sin2 in the basis.

The paired real/imaginary top-block energy equals Re(a* X a). Its oscillatory self-product uses w^2w^-2=1, cancelling apparent fourth harmonics. Polarization preserves this cancellation for every bilinear matrix entry. The mixed blocks involve f1, w^2f1, f2, w^-2f2 in the documented paired order. Their eta modes are only ±1, ±3 and xi modes only 0, ±2. Expansion gives the eight displayed columns with factor 1/sqrt(2). The lower mixed blocks must be transposes of the skew B blocks; both implementations do this.

The scalar space therefore has 30 X moments, 10 constant Z moments, four eta moments, and three xi moments: 47 independent functions. The 23 outputs select distinct coordinates. This gives 24 auxiliaries. The third reconstruction and full rank check additionally establish that the products span all 47 coordinates, not a smaller hidden subspace.

### Matching, indefinite matrices, and the axis

For delta>0, the displayed 2×2 coefficient system has A* A=delta I, and its right side has norm squared delta. Thus M is unitary and equals I at the base angle. Its first column uses only modes {0,-2}, and its second only {0,+2}. Multiplying the fixed Gram factor R0 preserves those supports. For delta=0, g=h=0 implies v=0; M=I is a separately valid constant choice with no division.

The V pairing is constant. The potentially cubic U term is multiplied by h·(w0/w)^2 and reduces to a first harmonic; the remaining row products also have only frequencies ±1. At the base angle the value agrees by Gram polarization. Differentiating the fixed Gram and constant V pairing removes the derivative of the moving factor. The U pairing has the required derivative. A real first harmonic is determined by its value and derivative at one angle.

This reasoning requires symmetry, not positivity, of S. Taking S=J_u D J_u^T is therefore valid for every symmetric D, including indefinite D and negative multiples of kernel projectors. The constant real map -D J_u^T preserves the separate column supports, so all the converse's matching tests actually lie in the restricted 36-dimensional space.

Independently of the supplied checker, `matching_direct_audit.py` expands the cleared-denominator pairing identity for six unrestricted real symmetric S entries, arbitrary tau/rho and eta/xi, and unit-circle Laurent variables for both angles. Its numerator vanishes identically, and it obtains the same separate supports. It also checks generic scaled unitarity. The zero-denominator branch is verified by the written argument and the replay's separate branch checks.

Applying an arbitrary feasible linear functional to those tests gives inequality (27). For fixed u,D,X,Z,y, the expression is a polynomial in z. Extending only this scalar inequality from nonaxis z to axis z is legitimate by continuity. It does not require continuous matching factors, representing measures, or closedness of a projection.

### Determinant identity and both cone implications

With d=det X, A=I4 tensor X, and hat J=I4 tensor adj X, one has A hat J=d I16. The four skew blocks give B^T hat J B=Phi_y(adj X). Right multiplication by the stated R gives the block-triangular matrix in equation (16). Both det R and det(d A) equal d^20. In the polynomial ring, d is nonzero (it equals 1 at I), so cancelling d^20 gives det L=p before any evaluation. Evaluating afterwards at singular X is valid. No factor is removed from p.

Since L(e)=I20 and the slice map is linear, p(te-x)=det(t I20-L(x)) exactly. This is monic of degree 20. All roots, their multiplicities, and zeros are retained. The spectral theorem then gives exact slice membership, including singular X.

The forward implication uses the explicit finite evaluation functional at a slice angle. Every auxiliary is a finite scalar or matrix multiple in equation (25); H=P^T L P is PSD. No diverging boundary witness is needed.

For the converse, positivity of constant top tests first gives X PSD. The argument uses every feasible functional, not merely evaluations. For v in ker X, taking D=t vv^T and z=y makes the quadratic X contribution vanish. The remaining inequality holds for both signs and all magnitudes of t, forcing (J_u v)^T Q(y)(J_u v)=0. The Gram identity makes this the sum of four real squares. Each square separately vanishes for every u,v. Consequently range(B_j) is contained in range(X) for all four j.

Taking D=X-dagger and z=y then gives Z >= Phi_y(X-dagger), since X-dagger X X-dagger=X-dagger. The already proved range relations give X X-dagger B_j=B_j. Expanding the four completed squares proves actual slice positivity. Neither the pseudoinverse inequality alone nor a closure operation is used. The argument also covers X=0 and the origin.

## Exact execution and adversarial checks

Python 3.14.7 and SymPy 1.14.0 were used in a task-local virtual environment. No float, numerical rank, tolerance, sampled feasibility, or numerical SDP solver supplied mathematical acceptance evidence.

| Check | Reproduced result | Evidence |
|---|---|---|
| Original four-stage replay | core, independent, matching, determinant all exit 0 | `tests/run_summary.json`, `evidence/sanitized/full-replay.log` |
| Original coefficient comparison | 47 matrices, 60912 full entries, 28512 directed contractions, 298 nonzero upper entries | `audits/independent_comparison.json` |
| Original matching checker | 42 named check groups, including both branches and singular tests | `tests/matching.log` |
| Original determinant checker | all 400 block-product entries, generic adjugate and Phi identities | `tests/determinant.log` |
| Third coefficient reconstruction | symbolic P^T L P from source definitions, no handoff imports; same full counts | `evidence/independent/third-coefficients-ranks.json` |
| Third bounded ranks/kernels | full rank47, auxiliary rank24; full/output common kernels zero; auxiliary common vector kernel span(e33) | same |
| Exact selected minors | auxiliary determinant1/16; full determinant sqrt(2)/32 | same |
| Direct generic-S matching | numerator identically zero; supports [-2,0] / [0,2]; scaled unitarity | `evidence/independent/direct-generic-matching.json` |
| New singular/root regressions | exact original-p/characteristic-polynomial identities of degree20 | `evidence/independent/singular-root-checks.json` |
| Corruption/optimization controls | all nine expected outcomes reproduced | `evidence/sanitized/historical/controls-summary.json` |
| Separate guard correction | optimized execution rejected at all four standalone entry points; normal full replay passes | `evidence/sanitized/hardening-summary.json` |

The third rank computation uses the full 666×47 upper-entry coefficient map, not just the supplied saved rank certificate. Common vector kernels are computed from sums M_k^T M_k; their kernels equal the intersections for real matrices. The auxiliary kernel e33 is not removable because the output Z22 coefficient acts on it. No unrestricted search or PSD-specific parameter elimination was performed.

The new feasible regression uses an off-diagonal rank-three X, its exact supported inverse, and a positive Schur residual. Its slice and noncardinal-angle finite witness have rank13. Its degree20 polynomial has seven zero roots, with factorization recorded in the JSON. The angle w=3/5+4i/5 exercises nonzero sin2 placements.

The pseudoinverse-only invalid point reproduces the negative principal minor -1/2. Direct polynomial expansion gives t^14(t-2)^2(t-1)^2(t^2-2), explicitly retaining the negative root -sqrt(2). X=0 with nonzero y is also rejected, while the origin gives t^20. These supplement the universal proof; they do not replace it.

## Verifier audit and separate correction

I read generator.py, verify.py, run_all.py, independent_coefficients.py, check_matching.py, and check_determinant_identity.py. The independent coefficient program reconstructs before reading primary JSON, compares complete coefficient dictionaries including omitted zeros, checks all matrix counts and coordinate selectors, and separately contracts the lower triangle. It does not import the primary generator. Neither program loads a saved PASS certificate to determine success. Exact minors and the identity combination support the rank/kernel conclusions before the report fields are assigned.

The matching checker verifies universal algebraic components and exact regressions; the determinant checker verifies the universal elimination product. Their concluding mathematical deductions still depend on the written proof. In particular, a printed `PASS` is not a formal theorem checker.

Deliberate corruptions of coefficient values, a missing matrix, source bytes, output selectors, a matching Gram identity, and a Phi sign were all rejected under ordinary execution. A mutated matching script prints PASS under `-O`, while the official runner refuses optimized execution. This is a real fail-open standalone interface outside documented usage, not evidence of a false theorem. `corrections/assertion-guards.patch` adds the existing runner's guard to all four standalone verification files. It changes no formula, coefficient, proof, or theorem. The guarded copy passes the normal replay and refuses all four optimized invocations.

