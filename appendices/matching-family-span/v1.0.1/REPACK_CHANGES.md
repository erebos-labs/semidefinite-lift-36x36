# Distribution correction 1.0.1

The prior v1.0.0 ZIP remains unchanged. This is a separate repack of that verified fixed witness, not a new version of the published 36 by 36 lift.

Changed existing files: README version/licensing/distribution explanation; PROOF version and licensing paragraph only; changelog; regenerated manifest. Mathematical sections 2 through 9 of PROOF are byte-identical, as are all source files, certificates, witness data, checker/search code and original verification reports.

Added: scoped root MIT license, exact pinned Apache-2.0 text, component scope, upstream attribution/notice, this change record, preservation provenance and SHA-256 checksums. No mathematical correction, new search, altered witness, removed assertion or public release is implied.

The original primary report identifies the historical Python environment. It is preserved, not relabeled as a fresh replay. The original search log and discovery program are retained history; the repack replay does not execute the search. Fresh extraction verification and privacy/provenance evidence are retained separately by the preparer and are not private raw logs embedded in this public distribution.

`MANIFEST.json` covers payload files except itself and `CHECKSUMS.sha256`, avoiding a checksum cycle. `CHECKSUMS.sha256` covers every file except itself, including the manifest. The external ZIP checksum identifies this exact distribution.
