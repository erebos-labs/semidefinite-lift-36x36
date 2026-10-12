# Exact 36 × 36 semidefinite lift

Public candidate **1.0.2** provides an explicit homogeneous real symmetric pencil over Q(sqrt(2)), with **23 outputs and 24 auxiliaries**, for the entire closed cone defined by the pinned degree-20 polynomial. Singular X and zero-root multiplicities are included.

The mathematical payload is unchanged from the audited v1.0.0 theorem and the hardened v1.0.1 handoff. This release changes privacy, distribution metadata, portable audit paths, and packaging. The matrices, coordinate ordering, source, verifier formulas, and assertion guards are preserved.

Read [THEOREM_PROOF.md](THEOREM_PROOF.md) for the complete theorem, explicit pencil, matching construction, both implications, and singular completion. [COORDINATES.md](COORDINATES.md) specifies every moment and coefficient encoding. [AUDIT.md](AUDIT.md) is a clearly labeled public abridgment of an independent model audit.

## Reproduce

Use Python 3.11 or 3.14; local reproduction uses Python 3.14.7. The tested dependency versions are SymPy 1.14.0 and mpmath 1.3.0.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt -c constraints.txt
.venv/bin/python run_all.py
.venv/bin/python public_controls.py
.venv/bin/python privacy_scan.py
```

`run_all.py` performs four exact stages: core coefficients/ranks/regressions, independent Laurent reconstruction, matching, and polynomial determinant elimination. `public_controls.py` checks direct entry points, optimized-execution refusal, deliberate corruption, and recomputation of forged saved certificates. It also reruns three supplemental reconstructions. The CI workflow runs these commands with read-only repository permissions; no hosted CI run is claimed for this local candidate.

Every revised verification entry point refuses `-O`, `-OO`, and positive PYTHONOPTIMIZE settings before claiming PASS. Verification requires assertions enabled. The removed historical scripts with the diagnosed hazard are retained privately, not included here.

All mathematical arithmetic is exact. There is no numerical feasibility solver, numerical rank threshold, or tolerance-based substitute for proof. The complete replay checks 47 matrices / 60,912 entries, 28,512 directed contractions, 298 nonzero upper entries, 42 matching check groups, the 400-entry elimination identity, and singular regressions. The bounded rank/kernel audit verifies ranks 47 and 24 and selected minor determinants sqrt(2)/32 and 1/16.

## Source and attribution

The source is OpenAI, *An Exact Semidefinite Lift of a Nonspectrahedral Hyperbolicity Cone* (5 October 2026), at commit [`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb). Unmodified proof sources, hashes and the Apache-2.0 license are under `source/`.

The source supplies the 100 × 100 / 307-auxiliary baseline, angular factorization, matching rotation, and singular-boundary converse. This construction restricts the matching columns to supports {0,-2} and {0,+2}, producing the smaller explicit pencil. It does not introduce a new general lifting principle.

**Nonspectrahedrality is inherited from the cited companion theorem.** This audit independently checked the exact lift, not that companion proof. See [ATTRIBUTION.md](ATTRIBUTION.md) and [NOTICE](NOTICE).

## Evidence and privacy

- `coefficients.json` and `generator.py`: generated matrices and exact implementation.
- `certificates/` and `audits/`: derived certificates, independent coefficient code, and mathematical review notes.
- `evidence/sanitized/`: explicitly labeled sanitized execution logs and public control summary. Filesystem roots are redacted; exit codes, commands and algebraic results are retained.
- `provenance/public-export.json`: source and mathematical-payload hashes, and the public edit scope.
- `SHA256SUMS.txt`: inventory for the packaged files.

Raw historical logs, private handoff archives, Library/account identifiers, operational notes, and personal filesystem paths remain private. Public logs are sanitized copies, never relabeled originals. To retain raw output from a new control run, set `SEMIDEFINITE_PRIVATE_RAW_DIR` to a directory outside the candidate. New runs may rewrite local result files.

The privacy scan is a deterministic pattern check supplemented by manual review; it is not a proof that every imaginable secret is absent. The release contains no binaries, account metadata, or historical unsafe verifier copies.

## Limits

This is model-assisted research with a written proof and exact finite verification. It is **not peer reviewed**, does not establish historical originality or priority, and makes no global minimality, Lean/formalization, numerical-conditioning, or solver-performance claim. The quantified theorem rests on the written argument; scripts do not formalize every logical deduction.

Human mathematical review can increase confidence but is not an unresolved proof obligation or a mandatory external-review gate identified by this audit. No specific mathematical gap was found.

## Licensing and release status

Original implementation code and new release materials use the [MIT license](LICENSE). Bundled upstream source and source-derived proof/review documentation retain Apache-2.0 terms. This is a scoped multi-license distribution: **MIT does not relicense third-party content**. See [LICENSING.md](LICENSING.md) for the file-level scope and ownership limitations.

This is a local candidate. No repository was created, no commit or push was made, and no publication or external message occurred. Proposed repository name: `semidefinite-lift-36x36`.
