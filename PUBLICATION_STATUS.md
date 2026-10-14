# Publication status and evidence scope

Status note prepared 15 October 2026. The main package was published as [v1.0.2](https://github.com/erebos-labs/semidefinite-lift-36x36/releases/tag/v1.0.2) at commit `2d2c170206ab937bcf56ff3290066cc3b2dfdd99`. The matching-family appendix was published as [matching-span-v1.0.1](https://github.com/erebos-labs/semidefinite-lift-36x36/releases/tag/matching-span-v1.0.1) at commit `9c48af3d19954d5ce797ee4b4c291e0f22bdea78`.

Local-candidate/no-publication wording inside the original reviewed packages describes their preparation stage. Immutable release assets retain that history; this status note identifies the later publication record. It changes no matrix, certificate, mathematical proof, or theorem scope.

The evidence consists of written mathematical arguments and exact computational checks. The main theorem includes singular cases. The appendix proves full span for the specified matching-function family, not a global lower bound on lift size. Nonspectrahedrality remains an inherited cited result whose proof was not independently audited in this lift review. No Lean verification, human peer review, global minimality, historical priority or solver-performance claim is made.

Publication receipts record hosted checks for the release commits. This status correction does not rerun those checks. A retrospective claim-level audit is in progress; neither this note nor a historical green badge constitutes a new verification under the internal publication standard adopted on 15 October 2026.

## Version-to-payload map

| Identifier | Meaning |
|---|---|
| theorem/generator 1.0.0 | retained mathematical payload headers |
| main distribution v1.0.2 | privacy/licensing/packaging release of that payload |
| matching-span-v1.0.1 | separate matching-family fixed-witness appendix |
| documentation revision 2026-10-15.1 | current landing-page status and source-inventory clarification; no theorem/code change |

Release objects and tags can change on hosting platforms. Use the exact commits above and source-asset hashes: main `71bcb0b6557f5a68e3970878088354fb7521ae780e8a7d402163a0efea6ca4e0`; appendix `54a6383455cc90329de0eb0851844890099f85340ba2f300cde4000131daf0b5`. The current source checksums describe this documentation revision; checksums inside the old released ZIPs remain unchanged.
