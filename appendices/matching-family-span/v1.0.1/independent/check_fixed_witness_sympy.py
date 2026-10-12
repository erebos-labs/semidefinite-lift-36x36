#!/usr/bin/env python3
"""Independently replay the single supplied span witness, using exact SymPy.

This performs no search and imports no other local generator or verifier.
It reconstructs actual Q(sqrt(2),i) matching coefficients from the source
formula, then checks the supplied rational scaling and exact determinant.
Dependency: SymPy (verified with 1.14.0).
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import sympy as s


HERE = Path(__file__).resolve().parent
PIN = "2d2c170206ab937bcf56ff3290066cc3b2dfdd99"
ROOT2 = s.sqrt(2)
I = s.I


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def exact(value):
    require(isinstance(value, (int, str)), "Only integer/rational-string inputs allowed")
    result = s.Rational(value)
    require(not result.atoms(s.Float), "Unexpected float")
    return result


def clean(value):
    return s.simplify(s.expand_complex(value))


def same(a, b, label):
    require(clean(a-b) == 0, label)


def zero_matrix(matrix, label):
    for entry in matrix:
        same(entry, 0, label)


def source_matching_column(record):
    u = s.Matrix([exact(x) for x in record["u"]])
    D = s.Matrix([[exact(x) for x in row] for row in record["D"]])
    require(u.shape == (4,1) and D.shape == (4,4), "Invalid source dimensions")
    require(D == D.T, "Matching matrix D must be exactly symmetric")
    c, sn, tau, rho = [exact(record[k]) for k in ("c", "s", "tau", "rho")]
    require(tau > 0, "The matching input must have tau>0")
    same(c*c+sn*sn, 1, "Angle pair must be on the unit circle")
    require(tau*c != 0 or tau*sn != 0, "Axis z is not admissible here")
    w0 = c+I*sn
    U0 = s.Matrix([[-1/w0, -w0], [sn/w0**2, sn*w0**2],
                    [c/w0**2, -c*w0**2]])/ROOT2
    V0 = s.Matrix([[0,0], [I,I], [1,1]])/ROOT2
    R0 = (tau*U0+rho*V0).applyfunc(clean)
    J = s.Matrix.hstack(s.Matrix(u[1:]), -u[0]*s.eye(3))
    S = J*D*J.T
    N = S*R0
    v = N[2,:]-I*N[1,:]
    g = clean(v[0]+s.conjugate(v[1]))
    h = clean(v[0]-s.conjugate(v[1]))
    delta = clean(g*s.conjugate(g)+h*s.conjugate(h))
    require(delta.is_Rational and delta >= 0, "Delta must be exact nonnegative rational")
    gh = [clean(ROOT2*g), clean(ROOT2*h)]
    for name, actual in zip(("g_hat","h_hat"), gh):
        saved = exact(record[name][0])+I*exact(record[name][1])
        same(actual, saved, "Incorrect recorded "+name)
    same(delta, exact(record["delta"]), "Incorrect source delta")
    same(2*delta, exact(record["delta_hat"]), "Incorrect normalized delta")

    K = -D*J.T
    pq = (K*R0).applyfunc(clean)
    p, q = pq[:,0], pq[:,1]
    if delta == 0:
        require(g == 0 and h == 0, "Zero delta branch has nonzero g or h")
        branch = "delta=0"
        A, B, C, Df = p, s.zeros(4,1), q, s.zeros(4,1)
    else:
        branch = "delta>0"
        r = clean(g*s.conjugate(g)/delta)
        t = clean(h*s.conjugate(h)/delta)
        vgh = clean(g*h/delta)
        phase = w0**2
        A = r*p-vgh*q
        B = phase*(t*p+vgh*q)
        C = s.conjugate(vgh)*p+r*q
        Df = s.conjugate(phase)*(-s.conjugate(vgh)*p+t*q)
    require(record["branch"] == branch, "Wrong denominator branch")
    actual_coordinates = []
    for vector in (A,B,C,Df):
        actual_coordinates.extend([clean(s.re(x)) for x in vector])
        actual_coordinates.extend([clean(s.im(x)) for x in vector])
    actual_coordinates.extend(list(u))
    require(len(actual_coordinates) == 36, "Wrong coordinate vector length")
    scaled = [clean(ROOT2*x) if row < 32 else x
              for row,x in enumerate(actual_coordinates)]
    require(all(x.is_Rational for x in scaled), "Scaled source column is not rational")
    return s.Matrix(scaled), branch


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--witness",type=Path,default=HERE/"span_witness.json")
    parser.add_argument("--output",type=Path,default=HERE/"sympy_witness_review.json")
    args=parser.parse_args()
    raw=args.witness.read_bytes()
    payload=json.loads(raw)
    require(payload["source_commit"] == PIN, "Wrong authoritative source pin")
    groups=["A.real","A.imag","B.real","B.imag","C.real","C.imag",
            "D.real","D.imag","u"]
    labels=[f"{group}.{i}" for group in groups for i in range(4)]
    require(payload["test_coordinate_labels"] == labels, "Wrong source coordinate order")
    records=payload["tests"]
    require(len(records) == len(payload["scaled_columns"]) == 36,
            "This verifier expects the single supplied 36-column witness")
    regenerated=[]
    counts=Counter()
    for j,record in enumerate(records):
        require(record["selected_column"] == j, "Column ordering mismatch")
        require(record["rank_after"] == j+1, "Inconsistent recorded prefix rank")
        if record["kind"] == "constant":
            k=record["index"]
            require(j < 4 and k == j, "Unexpected additional constant test")
            require(record["normalization"] == "1/sqrt(2)", "Constant scaling mismatch")
            column=s.zeros(36,1)
            column[k]=1  # Source vector e_k/sqrt(2), scaled by sqrt(2).
            counts["constant"]+=1
        else:
            require(record["kind"] == "matching", "Unknown test class")
            column,branch=source_matching_column(record)
            counts[branch]+=1
        saved=s.Matrix([exact(x) for x in payload["scaled_columns"][j]])
        require(saved.shape == (36,1), "Stored column dimension mismatch")
        zero_matrix(column-saved, f"Reconstructed column {j} differs from saved data")
        regenerated.append(column)
    matrix=s.Matrix.hstack(*regenerated)
    saved_rows=s.Matrix([[exact(x) for x in row] for row in payload["scaled_matrix_rows"]])
    zero_matrix(matrix-saved_rows, "Row-major and column-major matrices disagree")
    require(counts == Counter({"constant":4,"delta=0":8,"delta>0":24}),
            "Wrong composition of this fixed witness")
    require(payload["branches_selected"] == {"delta=0":8,"delta>0":24},
            "Wrong recorded branch counts")

    determinant=matrix.det(method="domain-ge")
    require(determinant == exact("-1024/125"), "Exact normalized determinant changed")
    same(determinant, exact(payload["scaled_matrix_determinant"]), "Wrong saved determinant")
    require(determinant != 0, "Witness is singular")
    require(payload["rank"] == 36, "Wrong claimed rank")
    original=determinant/(s.Integer(2)**16)
    same(original, exact("-1/8000"), "Incorrect source-coordinate determinant")
    same(original, exact(payload["original_matrix_determinant"]), "Wrong saved source determinant")

    result={
        "status":"PASS",
        "review_scope":"Only the supplied fixed 36-column witness; no new search",
        "source_commit":PIN,
        "witness_sha256":hashlib.sha256(raw).hexdigest(),
        "sympy_version":s.__version__,
        "arithmetic":"Exact source Q(sqrt(2),i) coefficients, rational normalized determinant",
        "implementation_independence":"Imports no local generator, search code, or other verifier",
        "tests_reconstructed":36,
        "column_entries_compared":36*36,
        "matching_tests":32,
        "constant_first_block_tests":4,
        "matching_branches":{"delta=0":8,"delta>0":24},
        "all_matching_inputs_nonaxis":True,
        "all_matching_D_exactly_symmetric":True,
        "coordinate_order_verified":labels,
        "scaled_determinant":str(determinant),
        "row_scaling_determinant":str(s.Integer(2)**16),
        "original_determinant":str(original),
        "rank":36,
        "conclusion":"The supplied augmented proof-test family spans the source 36-dimensional W",
        "limitations":["Does not prove global lift minimality or auxiliary minimality",
                       "Does not equate positivity on finitely many spanning vectors with PSD",
                       "Does not perform coefficient redundancy analysis"]}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
