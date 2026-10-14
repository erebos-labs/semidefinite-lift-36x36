# Matching-test span follow-up — v1.0.1

**Completed result:** the specified matching and constant tests span all 36
real test coordinates. An explicitly admissible 36-column witness has exact
determinant **−1/8000** in the published coefficient basis.

This is a separate follow-up to the user-designated published v1.0.2 at
erebos-labs/semidefinite-lift-36x36, pinned to
`2d2c170206ab937bcf56ff3290066cc3b2dfdd99`.
The existing 36 by 36 lift with 24 auxiliaries is preserved.

## Read and replay

Read `PROOF.md` for the complete theorem, explicit admissible parameter table,
coefficient-generation formulas, determinant and integer inverse certificates,
and the exact scope of the obstruction.

```bash
python verify_span.py
```

Python 3.10 or later and the standard library suffice. No network, solver,
floating-point arithmetic, or search is used. Run without `-O`.
The ordinary replay is read only; `--report path.json` writes a fresh report.

Expected core results:

```text
reconstructed_test_columns: 36
constant_tests: 4
matching_delta_zero: 8
matching_delta_positive: 24
integer_inverse_identities: W*U = U*W = 60*I36
integer_determinant: -491520
scaled_determinant: -1024/125
actual_determinant: -1/8000
rank: 36
status: PASS
```

## Files

| File | Purpose |
|---|---|
| `PROOF.md` | Complete proof and finite witness specification |
| `span_algebra.py` | Independent Gaussian-rational Laurent reconstruction |
| `verify_span.py` | Dependency-free fixed-witness replay |
| `exact_rank_certificate.json` | Integer witness, integer inverse multiplier, determinant and elimination data |
| `verification_report.json` | Completed primary verification record |
| `independent/span_witness.json` | All admissible inputs, branch values, and 1,296 rational matrix entries |
| `independent/check_fixed_witness_sympy.py` | Separate source-level exact verifier importing no local code |
| `independent/sympy_witness_review.json` | Completed independent verification record |
| `independent/declared_grid.json` | The original finite search bound |
| `independent/search_log.json` | Search stopped at candidate 638 of at most 1,881 |
| `independent/search_matching_span.py` | Retained discovery code; unnecessary for verification |
| `independent/README.md` | Search provenance and parameter table |
| `independent/MATCHING_FAMILY_AUDIT.md` | Independent formula and admissibility audit |
| `retrieval.json` | Exact source URLs, commit, byte sizes, and SHA-256 hashes |
| `source/` | Unchanged requested pinned files and their license |
| `MANIFEST.json` | SHA-256 inventory of the handoff files |

The optional separate replay requires SymPy, checked with version 1.14.0:

```bash
python independent/check_fixed_witness_sympy.py --output independent_review_new.json
```

The requested source files retain their original internal version headers.
Their pinned commit is authoritative; this follow-up does not edit them.
Original implementation and new distribution materials use scoped MIT terms
only to the extent their contributor holds the relevant rights. The bundled
source-derived theorem exposition and adapted mathematical notes retain
Apache-2.0 terms. `source/LICENSE` is the unchanged MIT file retrieved from the
baseline repository root; it does not relicense that exposition. The Apache
text is in `LICENSES/Apache-2.0.txt`. See `LICENSING.md`, `NOTICE` and
`ATTRIBUTION.md` for the component scope and upstream provenance.

This v1.0.1 repack corrects licensing and distribution metadata only. The
mathematics, fixed witness, rank/inverse certificates and checker code are
unchanged from v1.0.0. `MANIFEST.json` hashes all payload files except itself
and `CHECKSUMS.sha256`; `CHECKSUMS.sha256` also hashes the manifest and excludes
only itself. `REPACK_CHANGES.md` and `PUBLIC_PROVENANCE.json` record preservation.

## Scope and stopping

This proves that no proper **real linear test subspace** contains all of the
specific equation-(26) matching tests, with the published rotation rule and
both branches, together with the constant tests used for X positive
semidefinite. The finite list is a linear-span certificate; its 36 separate
quadratic inequalities are not asserted to replace a PSD condition.

No alternative lift, smaller sufficient certificate family, global minimality,
or novelty claim is investigated. No old coefficient-rank/common-kernel audit
is rerun. The lift, all 23 outputs and 24 auxiliaries, its degree-20 cone
definition, and its singular-boundary proof remain unchanged.

The following describes preparation of the original fixed-witness package, before its publication as `matching-span-v1.0.1`. See the repository-root [publication status](../../../PUBLICATION_STATUS.md) for the later publication record.

The search stopped upon the first exact full-rank witness. Subsequent work
only verified this fixed witness and prepared this handoff. No publication,
repository mutation, external communication, or scheduled follow-up was made.
