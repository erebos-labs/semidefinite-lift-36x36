# Component licensing scope

This distribution preserves the baseline project's scoped licensing. MIT applies only to original material to the extent the contributor owns the relevant original rights. It does not replace Apache-2.0 terms, claim ownership of mathematical facts, adjudicate authorship of machine-generated material, or grant rights in third-party content or trademarks.

| Component | Terms |
|---|---|
| `source/THEOREM_PROOF.md` | Apache-2.0, retaining the pinned baseline's terms for source-derived exposition. |
| `PROOF.md`, `independent/MATCHING_FAMILY_AUDIT.md`, and mathematical exposition in `independent/README.md` | Apache-2.0, conservatively retaining the source terms for adapted exposition and applying the same terms to any original protectable additions within the contributor's rights. |
| `span_algebra.py`, `verify_span.py`, `independent/check_fixed_witness_sympy.py`, `independent/search_matching_span.py`, `source/generator.py` | Scoped MIT for original implementation, subject to underlying third-party rights and notices. |
| `source/COORDINATES.md`, README, changelog, new distribution/attribution/licensing/provenance materials | Scoped MIT to the extent the contributor owns their original expression; any source-derived exposition retains applicable Apache-2.0 terms. |
| JSON witness, rank/inverse certificates, exact arithmetic data, recorded reports, source manifests and checksums | No exclusive rights are claimed in mathematical facts. Any protectable original organization contributed here is offered under scoped MIT; underlying third-party terms remain applicable. |
| `source/LICENSE`, root `LICENSE` | Unchanged MIT text from the pinned baseline root. It does not relicense the Apache components. |
| `LICENSES/Apache-2.0.txt` | Unmodified Apache-2.0 text inherited from the upstream pin. |
| Python, SymPy and other dependencies | Not vendored. Their own licenses apply. |

The upstream mathematical construction is OpenAI's *An Exact Semidefinite Lift of a Nonspectrahedral Hyperbolicity Cone* (5 October 2026), openai/math commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`. The directly consumed baseline is erebos-labs/semidefinite-lift-36x36 commit `2d2c170206ab937bcf56ff3290066cc3b2dfdd99`. The baseline's `LICENSING.md` and `NOTICE` establish the retained scope; their exact pin and hashes are recorded in `PUBLIC_PROVENANCE.json`.

Apache-2.0 text SHA-256: `c71d239df91726fc519c6eb72d318ec65820627232b2f796219e87dcf35d0ab4`. Preserve applicable licenses, attribution and modification notices on redistribution. This component map is a distribution statement, not a finding of legal ownership, patent clearance or protectability in every jurisdiction.
