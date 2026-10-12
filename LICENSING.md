# Licensing scope and provenance

This candidate uses MIT for original implementation code and new release materials, following the intended public project convention. That grant is limited to rights held by the project's contributor; it does not assert ownership of third-party content, mathematical facts, or copyright in purely machine-generated output.

## File-level terms

| Material | Terms |
|---|---|
| Original Python implementation/checkers, public release utilities, privacy/control/packaging scripts, CI configuration, new README and new release metadata | MIT, to the extent the contributor owns the relevant original rights. See root `LICENSE`. |
| Unmodified `source/build/`, `source/UPSTREAM_README.md`, and upstream `source/LICENSE` | Apache-2.0 from the pinned openai/math repository. No relicensing is claimed. |
| `THEOREM_PROOF.md`, source-derived mathematical review notes in `audits/*.md`, and `AUDIT.md` | Apache-2.0, conservatively retaining the source's terms for adapted exposition and licensing any original additions under the same terms. |
| Generated matrix/certificate data and source hash manifests | No claim of exclusive rights over mathematical facts. Any original protectable organization contributed here is available under MIT; underlying third-party notices remain applicable. |
| Dependency packages and CI actions | Referenced/installed, not vendored. Their own licenses apply. No rights in those projects are granted here. |

`source/LICENSE` is the unmodified Apache-2.0 text at commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`, Git blob `261eeb9e9f8b2b4b0d119366dda99c6fd7d35c64`, SHA-256 `c71d239df91726fc519c6eb72d318ec65820627232b2f796219e87dcf35d0ab4`. It was independently checked against that exact public revision. `LICENSES/Apache-2.0.txt` reproduces it for the adapted documentation.

The pinned repository root, the preprints ancestor, and the complete relevant source subtree contained no additional NOTICE or component license file. Source attribution is retained in the unchanged files and in `NOTICE` and `ATTRIBUTION.md`. The upstream files are not modified. The theorem note prominently identifies its compression/adaptation of the source in its existing introduction; `NOTICE` records that distinction. New release edits are documented in `PUBLIC_CHANGES.md`.

The redistribution basis is the pinned repository's Apache-2.0 grant and its redistribution conditions: include the license, preserve relevant notices/attribution, identify modifications, and retain any applicable upstream NOTICE. See the [official Apache-2.0 text, section 4](https://www.apache.org/licenses/LICENSE-2.0).

## Limits and recommendation

No conflicting component license or missing upstream NOTICE was found for these bundled source files. No license grant is made over cited papers, source trademarks, external services, or third-party content merely because it is cited.

Contributor ownership and the jurisdiction-specific copyright status of machine-generated additions were not legally adjudicated. The MIT convention identifies the intended original-code license; it does not certify authorship or cure third-party rights. Preserve the scoped terms and notices if this candidate is published. This is an evidence-backed packaging assessment, not a legal opinion about every possible patent, ownership, or jurisdictional issue.
