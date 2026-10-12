#!/usr/bin/env python3
"""Replay the fixed matching-test full-span certificate, without any search.

Usage: python verify_span.py [--report verification_report.json]
Python >= 3.10; standard library only.  No network, SDP solver, or floats.
The established pencil coefficient/kernel calculations are not executed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from math import prod
from pathlib import Path
import platform

from span_algebra import determinant_certificate, test_vector

PIN = "2d2c170206ab937bcf56ff3290066cc3b2dfdd99"
ROOT = Path(__file__).resolve().parent
LABELS = [f"{group}.{part}.{k}" for group in "ABCD"
          for part in ("real", "imag") for k in range(4)] + [f"u.{k}" for k in range(4)]


def read_json(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    if not __debug__:
        raise RuntimeError("Do not use -O: exact assertions must remain enabled")
    provenance = read_json(ROOT/"retrieval.json")
    assert provenance["commit"] == PIN
    for item in provenance["files"]:
        source = ROOT/"source"/item["path"]
        assert source.stat().st_size == item["bytes"], item["path"]
        assert digest(source) == item["sha256"], item["path"]
        assert f"/{PIN}/" in item["url"]
    witness_path = ROOT/"independent/span_witness.json"
    witness = read_json(witness_path)
    assert witness["source_commit"] == PIN
    assert witness["test_coordinate_labels"] == LABELS
    assert len(witness["tests"]) == 36
    assert [t.get("index") for t in witness["tests"][:4]] == list(range(4))
    columns, counts = [], Counter()
    for column_number, test in enumerate(witness["tests"]):
        if column_number < 4:
            assert test["kind"] == "constant"
            assert test["normalization"] == "1/sqrt(2)"
        else:
            assert test["kind"] == "matching"
        column, delta_hat, branch = test_vector(test)
        if delta_hat is not None:
            assert delta_hat == F(test["delta_hat"])
            assert delta_hat/2 == F(test["delta"])
            assert test["branch"] == ("delta>0" if branch == "positive" else "delta=0")
        counts[branch] += 1
        columns.append(column)
    assert dict(counts) == {"constant": 4, "zero": 8, "positive": 24}
    assert columns == [[F(x) for x in column] for column in witness["scaled_columns"]]
    matrix = [list(row) for row in zip(*columns)]
    assert matrix == [[F(x) for x in row] for row in witness["scaled_matrix_rows"]]
    # A separate integer fraction-free elimination checks the searcher's result.
    computed = determinant_certificate(matrix)
    saved = read_json(ROOT/"exact_rank_certificate.json")
    for key, value in computed.items():
        assert value == saved[key], f"Rank certificate mismatch: {key}"
    assert computed["scaled_determinant"] == "-1024/125"
    assert computed["actual_determinant"] == "-1/8000"
    assert witness["scaled_matrix_determinant"] == computed["scaled_determinant"]
    assert witness["original_matrix_determinant"] == computed["actual_determinant"]
    assert computed["integer_determinant"] == -491520
    assert prod(computed["column_denominators"]) == 60000
    # An integer matrix identity supplies a rank certificate independent of
    # trusting a determinant/rank routine. Check both multiplication orders.
    integer_matrix = computed["integer_matrix"]
    multiplier = saved["inverse_certificate"]["integer_multiplier"]
    identity_scale = saved["inverse_certificate"]["identity_scale"]
    assert identity_scale == 60
    assert len(multiplier) == 36 and all(len(row) == 36 for row in multiplier)
    assert all(type(x) is int for row in multiplier for x in row)
    for left, right in ((integer_matrix, multiplier), (multiplier, integer_matrix)):
        for i in range(36):
            for j in range(36):
                assert sum(left[i][k]*right[k][j] for k in range(36)) == 60*(i == j)
    # Check the saved bounded-search provenance, not additional candidates.
    search_log = read_json(ROOT/"independent/search_log.json")
    assert search_log["witness_sha256"] == digest(witness_path)
    assert search_log["declared_grid_sha256"] == digest(ROOT/"independent/declared_grid.json")
    assert search_log["status"] == "FULL_RANK_WITNESS"
    assert search_log["stopped_immediately_at_full_rank"] is True
    assert search_log["matching_candidates_scanned"] == 638
    assert search_log["maximum_matching_candidates"] == 1881
    return {
        "status": "PASS",
        "scope": "Fixed equation-(26) matching tests plus first-block constants",
        "source_commit": PIN,
        "source_files_hash_checked": len(provenance["files"]),
        "reconstructed_test_columns": 36,
        "constant_tests": counts["constant"],
        "matching_delta_zero": counts["zero"],
        "matching_delta_positive": counts["positive"],
        "source_realification_order_checked": True,
        "admissibility_checked_exactly": True,
        "rotation_unitarity_and_normalization_checked_exactly": True,
        "rotation_row_identities_checked_exactly": True,
        "Laurent_support_and_stored_coefficients_checked_exactly": True,
        "integer_inverse_identities": "W*U = U*W = 60*I36",
        "integer_determinant": computed["integer_determinant"],
        "scaled_determinant": computed["scaled_determinant"],
        "actual_determinant": computed["actual_determinant"],
        "rank": 36,
        "search_executed_by_this_verifier": False,
        "established_pencil_coefficient_or_kernel_checks_executed": False,
        "python_version": platform.python_version(),
        "dependencies": "Python standard library only",
        "witness_sha256": digest(witness_path),
        "rank_certificate_sha256": digest(ROOT/"exact_rank_certificate.json"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = verify()
    output = json.dumps(result, indent=2) + "\n"
    if args.report:
        args.report.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
