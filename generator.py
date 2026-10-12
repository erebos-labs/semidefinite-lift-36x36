#!/usr/bin/env python3
"""Exact 36 x 36 lift, v1.0.0. Python >=3.10; standard library only.

The 47 coordinates are REAL MOMENTS, not complex Fourier coefficients.
All entries live in Q(sqrt(2)). No floating point arithmetic is used.
See COORDINATES.md and THEOREM_PROOF.md for the mathematical specification.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import json
from pathlib import Path


@dataclass(frozen=True)
class Q2:
    """a + b sqrt(2), with canonical rational a,b."""
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "a", F(self.a))
        object.__setattr__(self, "b", F(self.b))

    @staticmethod
    def of(value):
        if isinstance(value, Q2):
            return value
        if isinstance(value, (int, F)):
            return Q2(value)
        raise TypeError("Only exact integers, Fractions, and Q2 are accepted")

    def __add__(self, other):
        other = Q2.of(other)
        return Q2(self.a + other.a, self.b + other.b)
    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-Q2.of(other))

    def __rsub__(self, other):
        return Q2.of(other) + (-self)

    def __mul__(self, other):
        other = Q2.of(other)
        return Q2(self.a*other.a + 2*self.b*other.b,
                  self.a*other.b + self.b*other.a)
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Q2.of(other)
        denominator = other.a**2 - 2*other.b**2
        if denominator == 0:
            raise ZeroDivisionError
        return self * Q2(other.a/denominator, -other.b/denominator)

    def __bool__(self):
        return bool(self.a or self.b)

    def encode(self):
        return [self.a.numerator, self.a.denominator,
                self.b.numerator, self.b.denominator]

    @staticmethod
    def decode(value):
        return Q2(F(value[0], value[1]), F(value[2], value[3]))


ZERO, ONE, INV_SQRT2 = Q2(), Q2(1), Q2(0, F(1, 2))
PAIRS = [(i, j) for i in range(4) for j in range(i, 4)]
MOMENT_LABELS = [f"X{i+1}{j+1}.{h}" for i, j in PAIRS
                 for h in ("1", "cos2", "sin2")]
MOMENT_LABELS += [f"Z{i+1}{j+1}.1" for i, j in PAIRS]
MOMENT_LABELS += ["eta.cos1", "eta.sin1", "eta.cos3", "eta.sin3",
                  "xi.1", "xi.cos2", "xi.sin2"]
OUTPUT_INDICES = list(range(0, 30, 3)) + list(range(30, 40)) + [40, 41, 44]
AUX_INDICES = [k for k in range(47) if k not in OUTPUT_INDICES]
TEST_GROUPS = ["A.real", "A.imag", "B.real", "B.imag",
               "C.real", "C.imag", "D.real", "D.imag", "u"]
TEST_LABELS = [f"{g}.{i}" for g in TEST_GROUPS for i in range(4)]


def zeros(rows, cols=None):
    return [[ZERO for _ in range(cols or rows)] for _ in range(rows)]


def B(r):
    """Source skew matrix B(r), with coordinates 0,1,2,3."""
    out = zeros(4)
    for k in range(3):
        out[0][k+1] = r[k]
        out[k+1][0] = -r[k]
    return out


def put_block(matrix, row, col, block, scalar=ONE):
    for i in range(4):
        for j in range(4):
            matrix[4*row+i][4*col+j] = scalar*block[i][j]


def pencil(moments):
    """Return H(Lambda) from the 47 ordered exact moment values."""
    if len(moments) != 47:
        raise ValueError("Expected 47 moment coordinates")
    m = [Q2.of(x) for x in moments]
    x0, xc, xs, z = [zeros(4) for _ in range(4)]
    for k, (i, j) in enumerate(PAIRS):
        for h, target in enumerate((x0, xc, xs)):
            target[i][j] = target[j][i] = m[3*k+h]
        z[i][j] = z[j][i] = m[30+k]
    h = zeros(36)
    for group in range(8):
        put_block(h, group, group, x0)
    put_block(h, 8, 8, z)
    # First complex pair: a = A + B exp(-2 i theta).
    # Second complex pair: b = C + D exp(+2 i theta).
    for row, col, block, sign in [
        (0, 2, xc, 1), (0, 3, xs, 1),
        (1, 2, xs, -1), (1, 3, xc, 1),
        (4, 6, xc, 1), (4, 7, xs, -1),
        (5, 6, xs, 1), (5, 7, xc, 1),
    ]:
        put_block(h, row, col, block, Q2(sign))
        put_block(h, col, row, block, Q2(sign))
    e1, e2, e3c, e3s, k0, kc, ks = m[40:47]
    ell = [
        [-e1, (e3s-e2)/2, (e3c+e1)/2+k0],
        [e2, (e3c-e1)/2+k0, -(e3s+e2)/2],
        [-e1, e2-ks, e1+kc],
        [-e2, kc, ks],
        [-e1, (e3s-e2)/2, -(e3c+e1)/2+k0],
        [-e2, (e1-e3c)/2+k0, -(e3s+e2)/2],
        [-e1, e2+ks, -e1+kc],
        [e2, kc, -ks],
    ]
    for group, r in enumerate(ell):
        cross = B([INV_SQRT2*x for x in r])
        put_block(h, group, 8, cross)
        # B(r) is skew; the global pencil is symmetric by transposition.
        transpose = [list(row) for row in zip(*cross)]
        put_block(h, 8, group, transpose)
    return h


def coefficient_matrices():
    result = []
    for k in range(47):
        moment = [ZERO]*47
        moment[k] = ONE
        result.append(pencil(moment))
    return result


def sparse(matrix):
    """Upper-triangular entries, zero-based indices, no implicit factor 2."""
    return [[i, j, value.encode()] for i, row in enumerate(matrix)
            for j, value in enumerate(row) if i <= j and value]


def coefficient_payload():
    mats = coefficient_matrices()
    return {
        "version": "1.0.0", "matrix_size": 36, "moment_dimension": 47,
        "output_dimension": 23, "auxiliary_dimension": 24,
        "field": "Q(sqrt(2))",
        "scalar_encoding": "[a_num,a_den,b_num,b_den] means a+b*sqrt(2)",
        "entry_encoding": "[row,col,scalar], upper triangular, zero based; reflect without scaling",
        "test_labels": TEST_LABELS, "moment_labels": MOMENT_LABELS,
        "output_indices": OUTPUT_INDICES, "auxiliary_indices": AUX_INDICES,
        "coefficients": [sparse(matrix) for matrix in mats],
    }


def output_pencil(outputs, auxiliary):
    """H = sum x_nu A_nu + sum v_j B_j, in the documented ordering."""
    if len(outputs) != 23 or len(auxiliary) != 24:
        raise ValueError("Expected 23 output and 24 auxiliary coordinates")
    moment = [ZERO]*47
    for indices, values in ((OUTPUT_INDICES, outputs), (AUX_INDICES, auxiliary)):
        for k, value in zip(indices, values):
            moment[k] = Q2.of(value)
    return pencil(moment)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("coefficients.json"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(coefficient_payload(), indent=2)+"\n")
    print(f"Wrote {args.output}: 47 exact symmetric 36x36 coefficient matrices")


if __name__ == "__main__":
    main()
