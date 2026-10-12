#!/usr/bin/env python3
"""Replay exact coefficient, rank/kernel, and singular-boundary checks.

No numerical tolerances, eigenvalue approximations, or SDP solver.
The mathematical implications are proved in THEOREM_PROOF.md. This script
checks finite algebraic data and explicit exact boundary examples.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

from generator import (Q2, ZERO, ONE, INV_SQRT2, PAIRS, MOMENT_LABELS,
                       OUTPUT_INDICES, AUX_INDICES, coefficient_matrices,
                       coefficient_payload, pencil, B, zeros)

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a, b):
    require(len(a[0]) == len(b), "matrix product dimensions")
    out = zeros(len(a), len(b[0]))
    for i, row in enumerate(a):
        for k, val in enumerate(row):
            if val:
                for j, v in enumerate(b[k]):
                    if v:
                        out[i][j] = out[i][j] + val*v
    return out


def identity(n):
    out = zeros(n)
    for i in range(n):
        out[i][i] = ONE
    return out


def eliminate(a):
    """Exact row echelon form, rank, pivots; determinant when square."""
    a = [row[:] for row in a]
    n, m = len(a), len(a[0])
    pivot_row, determinant, pivots = 0, ONE, []
    for col in range(m):
        k = next((r for r in range(pivot_row, n) if a[r][col]), None)
        if k is None:
            continue
        if k != pivot_row:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            determinant = -determinant
        v = a[pivot_row][col]
        determinant = determinant*v
        pivots.append([pivot_row, col, v.encode()])
        for r in range(pivot_row+1, n):
            if a[r][col]:
                ratio = a[r][col]/v
                for c in range(col, m):
                    a[r][c] = a[r][c] - ratio*a[pivot_row][c]
        pivot_row += 1
        if pivot_row == n:
            break
    if n == m and pivot_row != n:
        determinant = ZERO
    return {"rank": pivot_row, "determinant": determinant.encode() if n == m else None,
            "pivots": pivots}


def source_hashes():
    manifest = json.loads((ROOT/"source/manifest.json").read_text())
    require(manifest["commit"] == "fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb",
            "wrong source pin")
    for rec in manifest["files"]:
        actual = hashlib.sha256((ROOT/"source"/rec["path"]).read_bytes()).hexdigest()
        require(actual == rec["sha256"], "source hash mismatch: "+rec["path"])
    return {"commit": manifest["commit"], "files_checked": len(manifest["files"])}


def stored_coefficients():
    generated = coefficient_payload()
    saved = json.loads((ROOT/"coefficients.json").read_text())
    require(saved == generated, "stored coefficient artifact differs from generator")
    require(len(OUTPUT_INDICES) == 23 and len(AUX_INDICES) == 24, "coordinate counts")
    require(sorted(OUTPUT_INDICES+AUX_INDICES) == list(range(47)), "coordinate partition")
    mats = coefficient_matrices()
    for k, matrix in enumerate(mats):
        require(matrix == transpose(matrix), "nonsymmetric coefficient "+str(k))
    return {"matrices": 47, "size": 36,
            "nonzero_upper_entries": sum(len(x) for x in generated["coefficients"]),
            "coefficient_sha256": hashlib.sha256((ROOT/"coefficients.json").read_bytes()).hexdigest()}


def redundancy_audit():
    """The one bounded audit, of precisely this fixed 36/24 construction."""
    mats = coefficient_matrices()
    m = [ZERO]*47
    common_kernel_combination = []
    for k, (i, j) in enumerate(PAIRS):
        if i == j:
            m[3*k] = m[30+k] = ONE
            common_kernel_combination += [3*k, 30+k]
    require(pencil(m) == identity(36), "identity common-kernel certificate failed")

    # Each selected upper-triangle row isolates its corresponding auxiliary.
    # u's three nonzero spatial coordinates have indices 33,34,35.
    aux_rows = []
    for i, j in PAIRS:
        aux_rows += [(i, 8+j), (i, 12+j)]
    aux_rows += [(4, 34), (0, 34), (12, 34), (12, 35)]
    minor = [[mats[k][i][j] for k in AUX_INDICES] for i, j in aux_rows]
    expected_diag = [ONE]*20 + [INV_SQRT2/2]*2 + [INV_SQRT2]*2
    require(all(minor[i][j] == (expected_diag[i] if i == j else ZERO)
                for i in range(24) for j in range(24)), "auxiliary minor not diagonal")
    aux_result = eliminate(minor)
    require(aux_result["rank"] == 24, "auxiliary coefficient rank")
    require(aux_result["determinant"] == Q2(F(1,16)).encode(), "auxiliary determinant")

    # An additional certificate of faithfulness of the full moment map. This
    # uses the same matrices and audit, not an alternative-lift search.
    constant_rows = [(i,j) for i,j in PAIRS] + [(32+i,32+j) for i,j in PAIRS]
    full_rows = constant_rows+aux_rows+[(8,33),(12,33),(16,35)]
    full_minor = [[mats[k][i][j] for k in range(47)] for i,j in full_rows]
    full_result = eliminate(full_minor)
    require(full_result["rank"] == 47, "full coefficient map rank")

    # Common VECTOR kernel of the auxiliary subfamily: distinct from the
    # zero PARAMETER-map kernel just certified above.
    pc_indices = [1, 13, 22, 28]
    pc = [[sum((mats[k][i][j] for k in pc_indices), ZERO)
           for j in range(36)] for i in range(36)]
    target = zeros(36)
    for i in range(32):
        target[i][i] = ONE
    require(multiply(pc,pc) == target, "auxiliary top-space invertibility certificate")
    vector_minor_rows, vector_minor_cols = [11,12,8], [32,34,35]
    vector_minor = [[mats[45][i][j] for j in vector_minor_cols]
                    for i in vector_minor_rows]
    vector_det = eliminate(vector_minor)
    require(vector_det["determinant"] == (-INV_SQRT2/2).encode(),
            "auxiliary remaining-vector minor")
    require(all(not mats[k][i][33] for k in AUX_INDICES for i in range(36)),
            "auxiliary common vector e33")
    require([mats[34][i][33] for i in range(36)] ==
            [ONE if i==33 else ZERO for i in range(36)],
            "output Z acts on auxiliary-only vector kernel")
    return {
        "audit_scope": "Only common kernels and linear auxiliary redundancy of this pencil",
        "common_kernel_dimension": 0,
        "output_family_common_vector_kernel_dimension": 0,
        "identity_combination_moment_indices": common_kernel_combination,
        "identity_combination_weights": [ONE.encode()]*len(common_kernel_combination),
        "auxiliary_map_domain_dimension": 24, "auxiliary_map_rank": 24,
        "auxiliary_kernel_basis": [], "auxiliary_kernel_dimension": 0,
        "auxiliary_minor_rows": aux_rows, "auxiliary_minor_columns": AUX_INDICES,
        "auxiliary_minor_diagonal": [x.encode() for x in expected_diag],
        "auxiliary_minor_determinant": aux_result["determinant"],
        "auxiliary_family_common_vector_kernel_dimension": 1,
        "auxiliary_family_common_vector_kernel_basis_indices": [33],
        "auxiliary_family_common_vector_kernel_certificate": {
            "top_invertible_combination_indices": pc_indices,
            "top_combination_squared": "diag(I32,0_4)",
            "remaining_coefficient_index": 45,
            "remaining_minor_rows": vector_minor_rows,
            "remaining_minor_columns": vector_minor_cols,
            "remaining_minor_determinant": vector_det["determinant"],
            "output_coefficient_acting_nontrivially_on_e33": 34,
            "removable_from_full_pencil": False,
        },
        "full_map_rank": 47, "full_map_kernel_dimension": 0,
        "full_minor_rows": full_rows, "full_minor_columns": list(range(47)),
        "full_minor_determinant": full_result["determinant"],
        "full_minor_elimination_pivots": full_result["pivots"],
        "removed_matrix_directions": 0, "removed_auxiliary_directions": 0,
        "global_minimality_claim": False,
    }


def theta_zero_test_matrix():
    """20x36 P(0) from the explicit real test basis."""
    p = zeros(20,36)
    for i in range(4):
        for block, groups in enumerate(((0,2),(1,3),(4,6),(5,7))):
            for group in groups:
                p[4*block+i][4*group+i] = ONE
        p[16+i][32+i] = ONE
    return p


def axis_slice(x_diagonal, z_diagonal):
    """Exact L_0(X,Z,eta=0,xi=1), Q(0,0,1)=diag(0,1,1)."""
    r = [[ZERO,ZERO,INV_SQRT2], [ZERO,INV_SQRT2,ZERO]]*2
    b = [B(col) for col in r]
    l = zeros(20)
    for j in range(4):
        for i in range(4):
            l[4*j+i][4*j+i] = Q2(x_diagonal[i])
            for k in range(4):
                l[4*j+i][16+k] = b[j][i][k]
                l[16+k][4*j+i] = b[j][i][k]
    for i in range(4):
        l[16+i][16+i] = Q2(z_diagonal[i])
    return l, b


def singular_tests():
    # A nonzero-y singular boundary point, with zero Schur residual.
    xd, zd = [1,0,1,1], [2,0,1,1]
    l, b = axis_slice(xd, zd)
    gram = zeros(16,20)
    for j in range(4):
        for i in range(4):
            gram[4*j+i][4*j+i] = Q2(xd[i])
            for k in range(4):
                gram[4*j+i][16+k] = b[j][i][k]
    require(multiply(transpose(gram), gram) == l, "singular slice exact Gram")
    require(eliminate(l)["rank"] == 12, "singular slice rank")
    m = [ZERO]*47
    for k,(i,j) in enumerate(PAIRS):
        if i == j:
            m[3*k] = m[3*k+1] = Q2(xd[i])  # theta=0
            m[30+k] = Q2(zd[i])
    m[44] = m[45] = ONE
    h = pencil(m)
    p = theta_zero_test_matrix()
    require(h == multiply(transpose(p), multiply(l,p)), "angular evaluation H=P^T L P")
    g = multiply(gram,p)
    require(h == multiply(transpose(g),g), "singular lift finite exact Gram")
    require(eliminate(h)["rank"] == 12, "singular lift rank")
    # Pseudoinverse inequality alone must NOT count as cone membership.
    # X=diag(1,0,0,0), Z=Phi_y(Xdag)=diag(0,0,1,1) has a failed range condition.
    bad, _ = axis_slice([1,0,0,0], [0,0,1,1])
    neg_minor = bad[3][3]*bad[16][16] - bad[3][16]*bad[16][3]
    require(neg_minor == Q2(-F(1,2)), "failed-range negative principal minor")
    require(pencil([ZERO]*47) == zeros(36), "origin witness")
    return {
        "singular_point": {"X_diagonal": xd,"Z_diagonal": zd,"y":[0,0,1]},
        "angular_witness": {"theta":0,"eta":0,"xi":1},
        "slice_rank":12,"lift_rank":12,"exact_gram_factorization":True,
        "all_four_range_conditions_checked_by_factorization":True,
        "pseudoinverse_only_counterexample": {
            "X_diagonal":[1,0,0,0],"Z_diagonal":[0,0,1,1],"y":[0,0,1],
            "negative_slice_minor_indices":[3,16],
            "negative_slice_minor_determinant":neg_minor.encode()},
        "origin_exact_zero_witness":True,
        "role":"Regression examples; universal boundary proof is in THEOREM_PROOF.md",
    }


def independent_coefficient_check():
    path = ROOT/"audits/independent_coefficients.py"
    require(path.exists(), "missing independently written Laurent reconstruction")
    import importlib.util
    spec = importlib.util.spec_from_file_location("independent_coefficients",path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    # The independent reconstruction exports its own complete coefficient
    # payload; the adapter below is completed against that documented format.
    if not hasattr(module, "compare_with_primary"):
        raise AssertionError("Independent comparator interface not yet supplied")
    return module.compare_with_primary(ROOT/"coefficients.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-independent",action="store_true",
                        help="Partial diagnostic run; does not establish full verification")
    parser.add_argument("--output-dir",type=Path,default=ROOT/"certificates")
    args = parser.parse_args()
    started=time.monotonic()
    result={"version":"1.0.0","arithmetic":"exact Q(sqrt(2)), no floats",
            "python":platform.python_version(),"checks":{}}
    checks=[("source_hashes",source_hashes),
            ("stored_coefficients",stored_coefficients),
            ("bounded_redundancy_audit",redundancy_audit),
            ("singular_boundary_regressions",singular_tests)]
    if not args.skip_independent:
        checks.append(("independent_laurent_reconstruction",independent_coefficient_check))
    for name,fn in checks:
        result["checks"][name]=fn()
        print("PASS",name,flush=True)
    result["coefficient_verification_complete"]=not args.skip_independent
    result["universal_matching_checker_included_here"]=False
    result["universal_matching_checker_command"]="python audits/check_matching.py"
    result["elapsed_seconds_for_information_only"]=round(time.monotonic()-started,3)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/"verification.json").write_text(json.dumps(result,indent=2)+"\n")
    (args.output_dir/"rank_kernel.json").write_text(json.dumps(result["checks"]["bounded_redundancy_audit"],indent=2)+"\n")
    print("Exact algebraic checks completed"+(" (PARTIAL: independent comparison skipped)" if args.skip_independent else ""))


if __name__ == "__main__":
    main()
