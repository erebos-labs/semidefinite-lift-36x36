#!/usr/bin/env python3
"""Independent exact expansion of q_i^T L_theta q_j and one redundancy audit.

Only Python's standard library is used.  The pencil is built directly from
the source Laurent formulas U,V and the source realification order, not from
the compact block formula in the primary generator.  Every arithmetic
operation takes place in Q(sqrt(2),i)[w,w^-1].  The final matrices are in
Sym^36(Q(sqrt(2))).  No numerical tolerance or sampling is used.

Run from any directory:
    python3 independent_coefficients.py --out-dir PATH

JSON stores each real coefficient a+b*sqrt(2) as the two rational strings
[a,b].  Matrix entries are sparse upper-triangular entries; lower entries
are their symmetric copies.  All indices in generated certificates are
zero based.  There are 47 coefficient matrices in the requested order.
"""

from __future__ import annotations

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path


@dataclass(frozen=True)
class Q2:
    a: F = F(0)
    b: F = F(0)

    def __add__(self, other):
        if not isinstance(other, Q2):
            other = Q2(F(other))
        return Q2(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-other if isinstance(other, Q2) else -Q2(F(other)))

    def __mul__(self, other):
        if not isinstance(other, Q2):
            other = Q2(F(other))
        return Q2(self.a * other.a + 2 * self.b * other.b,
                  self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, Q2):
            other = Q2(F(other))
        denominator = other.a * other.a - 2 * other.b * other.b
        if denominator == 0:
            raise ZeroDivisionError
        return self * Q2(other.a / denominator, -other.b / denominator)

    def __bool__(self):
        return bool(self.a or self.b)

    def json(self):
        return [str(self.a), str(self.b)]


ZERO = Q2()
ONE = Q2(F(1))
SQRT2_OVER2 = Q2(F(0), F(1, 2))


@dataclass(frozen=True)
class C2:
    """Complexification of Q(sqrt(2))."""
    re: Q2 = ZERO
    im: Q2 = ZERO

    def __add__(self, other):
        return C2(self.re + other.re, self.im + other.im)

    def __neg__(self):
        return C2(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        return C2(self.re * other.re - self.im * other.im,
                  self.re * other.im + self.im * other.re)

    def conjugate(self):
        return C2(self.re, -self.im)

    def __bool__(self):
        return bool(self.re or self.im)


CZ = C2()
CO = C2(ONE)
CI = C2(ZERO, ONE)


def clean(p):
    return {k: v for k, v in p.items() if v}


def lp_const(x):
    if not isinstance(x, C2):
        x = C2(x if isinstance(x, Q2) else Q2(F(x)))
    return {0: x} if x else {}


def lp_add(*polys):
    result = {}
    for p in polys:
        for k, v in p.items():
            result[k] = result.get(k, CZ) + v
    return clean(result)


def lp_scale(p, x):
    if not isinstance(x, C2):
        x = C2(x if isinstance(x, Q2) else Q2(F(x)))
    return clean({k: v * x for k, v in p.items()})


def lp_mul(p, q):
    result = {}
    for k, v in p.items():
        for ell, value in q.items():
            result[k + ell] = result.get(k + ell, CZ) + v * value
    return clean(result)


def lp_shift(p, n):
    return {k + n: v for k, v in p.items()}


def lp_conjugate(p):
    return {-k: v.conjugate() for k, v in p.items()}


def lp_real(p):
    return lp_scale(lp_add(p, lp_conjugate(p)), F(1, 2))


def lp_imag(p):
    return lp_scale(lp_add(p, lp_scale(lp_conjugate(p), -1)),
                    C2(ZERO, Q2(F(-1, 2))))


def cosine(k):
    return {k: C2(Q2(F(1, 2))), -k: C2(Q2(F(1, 2)))}


def sine(k):
    return {k: C2(ZERO, Q2(F(-1, 2))),
            -k: C2(ZERO, Q2(F(1, 2)))}


UPPER = [(i, j) for i in range(4) for j in range(i, 4)]
TEST_GROUPS = ['A_re', 'A_im', 'B_re', 'B_im',
               'C_re', 'C_im', 'D_re', 'D_im', 'u']


def source_pencil_parts():
    """Sparse 20x20 coefficient matrices, directly from angular/slices.tex."""
    parts = {}
    for name in ['X', 'Z']:
        for i, j in UPPER:
            matrix = {}
            offsets = [0, 4, 8, 12] if name == 'X' else [16]
            for offset in offsets:
                matrix[offset+i, offset+j] = lp_const(1)
                matrix[offset+j, offset+i] = lp_const(1)
            parts[name, i, j] = matrix

    c, s = cosine(1), sine(1)
    # U = 1/sqrt(2) * [[-w^-1,-w],
    #                       [s w^-2,s w^2], [c w^-2,-c w^2]].
    U = [[{-1: -CO}, {1: -CO}],
         [lp_shift(s, -2), lp_shift(s, 2)],
         [lp_shift(c, -2), lp_scale(lp_shift(c, 2), -1)]]
    V = [[{}, {}], [lp_const(CI), lp_const(CI)],
         [lp_const(1), lp_const(1)]]
    for scalar, complex_factor in [('eta', U), ('xi', V)]:
        factor = [[lp_scale(p, SQRT2_OVER2) for p in row]
                  for row in complex_factor]
        # Source order: Re first, Im first, Re second, Im second.
        real_columns = []
        for column in range(2):
            real_columns.append([lp_real(factor[k][column]) for k in range(3)])
            real_columns.append([lp_imag(factor[k][column]) for k in range(3)])
        matrix = {}
        for block, vector in enumerate(real_columns):
            # B(r) = [[0, r^T], [-r, 0]].
            for component, poly in enumerate(vector, start=1):
                if poly:
                    matrix[4*block, 16+component] = poly
                    matrix[16+component, 4*block] = poly
                    matrix[4*block+component, 16] = lp_scale(poly, -1)
                    matrix[16, 4*block+component] = lp_scale(poly, -1)
        parts[scalar] = matrix
    return parts


def test_basis():
    """Expand coefficients of a=A+B*w^-2, b=C+D*w^2 in source order."""
    real_basis = []
    for group in ['A', 'B', 'C', 'D']:
        shift = {'A': 0, 'B': -2, 'C': 0, 'D': 2}[group]
        offset = 0 if group in ['A', 'B'] else 8
        for scalar in [CO, CI]:
            monomial = {shift: scalar}
            real, imag = lp_real(monomial), lp_imag(monomial)
            for coordinate in range(4):
                components = {}
                if real:
                    components[offset + coordinate] = real
                if imag:
                    components[offset + 4 + coordinate] = imag
                real_basis.append(components)
    for coordinate in range(4):
        real_basis.append({16+coordinate: lp_const(1)})
    assert len(real_basis) == 36
    return real_basis


def moment_order():
    moments = []
    for i, j in UPPER:
        for trig in ['1', 'cos2', 'sin2']:
            moments.append(('X', i, j, trig))
    for i, j in UPPER:
        moments.append(('Z', i, j, '1'))
    moments.extend(('eta', trig) for trig in ['cos1', 'sin1', 'cos3', 'sin3'])
    moments.extend(('xi', trig) for trig in ['1', 'cos2', 'sin2'])
    assert len(moments) == 47
    return moments


def real_trig_coefficients(p):
    """Laurent-to-real-trigonometric conversion with exact reality checks."""
    assert lp_conjugate(p) == p, 'q_i^T L q_j must be pointwise real'
    result = {}
    if 0 in p:
        assert not p[0].im
        result['1'] = p[0].re
    for k in sorted(k for k in p if k > 0):
        # Coefficient of w^k is (cos coefficient - i sin coefficient)/2.
        if p[k].re:
            result['cos'+str(k)] = p[k].re * 2
        if p[k].im:
            result['sin'+str(k)] = p[k].im * (-2)
    return result


def generate():
    parts, basis, moments = source_pencil_parts(), test_basis(), moment_order()
    lookup = {m: k for k, m in enumerate(moments)}
    matrices = [{} for _ in moments]
    observed = {str(variable): set() for variable in parts}
    entry_expansions = 0
    for i in range(36):
        for j in range(i, 36):
            for variable, part in parts.items():
                products = []
                for row, qi in basis[i].items():
                    for column, qj in basis[j].items():
                        if (row, column) in part:
                            products.append(lp_mul(lp_mul(qi, part[row, column]), qj))
                expression = lp_add(*products)
                observed[str(variable)].update(expression)
                for trig, coefficient in real_trig_coefficients(expression).items():
                    coordinate = (variable, trig) if isinstance(variable, str) else (*variable, trig)
                    assert coordinate in lookup, ('unexpected frequency', coordinate, i, j)
                    matrices[lookup[coordinate]][i, j] = coefficient
                entry_expansions += 1
    # Verify the lower triangle by a second directed contraction.
    for i in range(36):
        for j in range(i):
            for variable, part in parts.items():
                products = []
                for row, qi in basis[i].items():
                    for column, qj in basis[j].items():
                        if (row, column) in part:
                            products.append(lp_mul(lp_mul(qi, part[row, column]), qj))
                directed = real_trig_coefficients(lp_add(*products))
                prefix = (variable,) if isinstance(variable, str) else variable
                expected_directed = {
                    moment[-1]: matrices[k][j, i]
                    for k, moment in enumerate(moments)
                    if moment[:-1] == prefix and (j, i) in matrices[k]
                }
                # Compare complete dictionaries, rejecting extra modes in
                # addition to checking every expected coefficient.
                assert directed == expected_directed, (
                    'lower-triangle mismatch or unexpected mode', variable, i, j,
                    directed, expected_directed)
                entry_expansions += 1
    expected = {}
    for variable in parts:
        if isinstance(variable, tuple):
            expected[str(variable)] = {-2, 0, 2} if variable[0] == 'X' else {0}
        else:
            expected[str(variable)] = {-3, -1, 1, 3} if variable == 'eta' else {-2, 0, 2}
    assert observed == expected, (observed, expected)
    return moments, matrices, {k: sorted(v) for k, v in observed.items()}, entry_expansions


def determinant(matrix):
    """Exact Gaussian elimination in Q(sqrt(2)); no numerical pivot tests."""
    a = [row[:] for row in matrix]
    n, result = len(a), ONE
    assert all(len(row) == n for row in a)
    for column in range(n):
        pivot = next((r for r in range(column, n) if a[r][column]), None)
        if pivot is None:
            return ZERO
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            result = -result
        value = a[column][column]
        result = result * value
        for row in range(column+1, n):
            if a[row][column]:
                multiplier = a[row][column] / value
                for k in range(column, n):
                    a[row][k] = a[row][k] - multiplier * a[column][k]
    return result


def extract_minor(matrices, row_positions, columns):
    return [[matrices[column].get(tuple(sorted(position)), ZERO)
             for column in columns] for position in row_positions]


def redundancy_audit(moments, matrices):
    """Exactly one bounded audit of this fixed pencil, with explicit minors."""
    output_columns = [3*k for k in range(10)] + list(range(30, 40)) + [40, 41, 44]
    auxiliary_columns = [k for k in range(47) if k not in output_columns]
    assert len(output_columns) == 23 and len(auxiliary_columns) == 24
    assert sorted(output_columns + auxiliary_columns) == list(range(47))

    # The sum of constant X/Z diagonal coefficients is exactly I_36.
    identity_columns = [0, 12, 21, 27, 30, 34, 37, 39]
    for i in range(36):
        for j in range(i, 36):
            value = sum((matrices[k].get((i, j), ZERO) for k in identity_columns), ZERO)
            assert value == (ONE if i == j else ZERO)

    auxiliary_rows = []
    for i, j in UPPER:
        auxiliary_rows.extend([(i, 8+j), (i, 12+j)])
    auxiliary_rows.extend([(4, 34), (0, 34), (12, 34), (12, 35)])
    auxiliary_minor = extract_minor(matrices, auxiliary_rows, auxiliary_columns)
    expected_diagonal = [ONE]*20 + [SQRT2_OVER2/2]*2 + [SQRT2_OVER2]*2
    assert all(auxiliary_minor[i][j] == (expected_diagonal[i] if i == j else ZERO)
               for i in range(24) for j in range(24))
    auxiliary_determinant = determinant(auxiliary_minor)
    assert auxiliary_determinant == Q2(F(1, 16))

    full_rows = []
    for i, j in UPPER:
        full_rows.extend([(i, j), (i, 8+j), (i, 12+j)])
    full_rows.extend((32+i, 32+j) for i, j in UPPER)
    full_rows.extend([(8, 33), (12, 33), (4, 34), (0, 34),
                      (20, 34), (12, 34), (12, 35)])
    full_minor = extract_minor(matrices, full_rows, list(range(47)))
    full_determinant = determinant(full_minor)
    assert full_determinant == Q2(F(0), F(1, 32))

    # Distinguish common kernels of coefficient MATRICES (vectors in R^36)
    # from the parameter-map kernel (vectors in R^24).  The auxiliary
    # coefficient family has common vector kernel span(e_33), even though
    # the full/output families have zero common kernel and the auxiliary
    # PARAMETER map is injective.
    auxiliary_top_identity_columns = [1, 13, 22, 28]
    def coefficient_entry(k, i, j):
        return matrices[k].get(tuple(sorted((i, j))), ZERO)
    top_annihilator = [
        [sum((coefficient_entry(k, i, j) for k in auxiliary_top_identity_columns), ZERO)
         for j in range(36)]
        for i in range(36)
    ]
    # This sum pairs A with B and C with D, and squares to diag(I_32,0_4).
    for i in range(36):
        for j in range(36):
            square_entry = sum((top_annihilator[i][k]*top_annihilator[k][j]
                                for k in range(36)), ZERO)
            assert square_entry == (ONE if i == j and i < 32 else ZERO)
    # Thus a common-kernel vector has zero top 32 coordinates.  The xi.cos2
    # coefficient then constrains u0,u2,u3 via this nonsingular minor.
    auxiliary_bottom_rows = [11, 12, 8]
    auxiliary_bottom_columns = [32, 34, 35]
    auxiliary_bottom_minor = [
        [coefficient_entry(45, i, j) for j in auxiliary_bottom_columns]
        for i in auxiliary_bottom_rows
    ]
    auxiliary_bottom_expected = [-SQRT2_OVER2, SQRT2_OVER2, SQRT2_OVER2]
    assert all(auxiliary_bottom_minor[i][j] ==
               (auxiliary_bottom_expected[i] if i == j else ZERO)
               for i in range(3) for j in range(3))
    auxiliary_bottom_determinant = determinant(auxiliary_bottom_minor)
    assert auxiliary_bottom_determinant == Q2(F(0), F(-1, 4))
    # The only surviving direction e_33 is actually killed by every aux
    # coefficient, proving equality of the claimed common kernel.
    assert all(not coefficient_entry(k, i, 33)
               for k in auxiliary_columns for i in range(36))
    # The constant Z22 OUTPUT coefficient acts nontrivially on e_33; the
    # auxiliary common vector kernel is therefore not a removable pencil row.
    assert coefficient_entry(34, 33, 33) == ONE

    return {
        'audit_count': 1,
        'pencil_dimension': 36,
        'moment_dimension': 47,
        'output_dimension': 23,
        'auxiliary_dimension': 24,
        'output_columns': output_columns,
        'auxiliary_columns': auxiliary_columns,
        'output_selector_rank': 23,
        'output_selector_certificate': 'Its output_columns minor is I_23.',
        'common_kernel_dimension': 0,
        'common_kernel_basis': [],
        'full_coefficient_family_common_vector_kernel_dimension': 0,
        'output_coefficient_family_common_vector_kernel_dimension': 0,
        'identity_combination_columns': identity_columns,
        'identity_combination_weights': [1]*8,
        'identity_combination_value': 'I_36',
        'auxiliary_map_rank': 24,
        'auxiliary_map_kernel_dimension': 0,
        'auxiliary_map_kernel_basis': [],
        'auxiliary_coefficient_family_common_vector_kernel_dimension': 1,
        'auxiliary_coefficient_family_common_vector_kernel_basis': [
            {'coordinate': 33, 'value': ['1', '0'], 'label': 'u[1]'}
        ],
        'auxiliary_top_annihilator_columns': auxiliary_top_identity_columns,
        'auxiliary_top_annihilator_weights': [1]*4,
        'auxiliary_top_annihilator_square': 'diag(I_32,0_4)',
        'auxiliary_bottom_constraint_coefficient': 45,
        'auxiliary_bottom_constraint_minor_rows': auxiliary_bottom_rows,
        'auxiliary_bottom_constraint_minor_columns': auxiliary_bottom_columns,
        'auxiliary_bottom_constraint_minor_diagonal': [q.json() for q in auxiliary_bottom_expected],
        'auxiliary_bottom_constraint_minor_determinant': auxiliary_bottom_determinant.json(),
        'auxiliary_common_vector_kernel_not_removable_reason':
            'The output coefficient H_34 (constant Z22) maps e_33 to e_33.',
        'auxiliary_minor_rows': auxiliary_rows,
        'auxiliary_minor_columns': auxiliary_columns,
        'auxiliary_minor_diagonal': [q.json() for q in expected_diagonal],
        'auxiliary_minor_determinant': auxiliary_determinant.json(),
        'full_coefficient_map_rank': 47,
        'full_coefficient_map_kernel_dimension': 0,
        'full_minor_rows': full_rows,
        'full_minor_columns': list(range(47)),
        'full_minor_determinant': full_determinant.json(),
        'scope': 'Only the coefficient common kernel and the linear auxiliary coefficient map with outputs fixed. No unrestricted search or minimality claim.'
    }


def json_dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')


def compare_with_primary(path):
    """Reconstruct first, then compare every coefficient with primary JSON.

    The primary generator is never imported or called.  Its JSON is read
    only after this independent source-level Laurent reconstruction.
    """
    moments, independent, support, contractions = generate()
    path = Path(path)
    primary = json.loads(path.read_text())
    assert primary['matrix_size'] == 36
    assert primary['moment_dimension'] == 47
    assert primary['output_dimension'] == 23
    assert primary['auxiliary_dimension'] == 24
    assert len(primary['coefficients']) == 47
    expected_test_labels = []
    for group in TEST_GROUPS:
        if group == 'u':
            expected_test_labels.extend('u.'+str(i) for i in range(4))
        else:
            letter, component = group.split('_')
            component = 'real' if component == 're' else 'imag'
            expected_test_labels.extend(f'{letter}.{component}.{i}' for i in range(4))
    expected_moment_labels = []
    for moment in moments:
        if moment[0] in ['X', 'Z']:
            name, i, j, trig = moment
            expected_moment_labels.append(f'{name}{i+1}{j+1}.{trig}')
        else:
            expected_moment_labels.append('.'.join(moment))
    assert primary['test_labels'] == expected_test_labels
    assert primary['moment_labels'] == expected_moment_labels
    expected_outputs = [3*k for k in range(10)] + list(range(30, 40)) + [40, 41, 44]
    assert primary['output_indices'] == expected_outputs
    assert primary['auxiliary_indices'] == [k for k in range(47) if k not in expected_outputs]
    for k, entries in enumerate(primary['coefficients']):
        reconstructed_primary = {}
        for i, j, value in entries:
            assert 0 <= i <= j < 36
            assert (i, j) not in reconstructed_primary, 'duplicate primary matrix entry'
            assert len(value) == 4
            coefficient = Q2(F(value[0], value[1]), F(value[2], value[3]))
            assert coefficient, 'sparse entries should be nonzero'
            reconstructed_primary[i, j] = coefficient
        # Dictionary equality checks both absent zeros and every nonzero.
        assert independent[k] == reconstructed_primary, ('coefficient mismatch', k, moments[k])
    return {
        'status': 'PASS',
        'method': 'Independent Laurent q_i^T L_theta q_j reconstruction compared with compact primary coefficient JSON',
        'arithmetic': 'Exact Q(sqrt(2),i); no floating point',
        'coefficient_matrices_compared': 47,
        'matrix_size': 36,
        'full_matrix_scalar_entries_compared': 47*36*36,
        'upper_triangle_scalar_entries_compared': 47*36*37//2,
        'nonzero_upper_entries': sum(len(matrix) for matrix in independent),
        'directed_source_contractions_checked': contractions,
        'frequency_support_checked': True,
        'lower_triangle_complete_coefficient_dictionary_equality': True,
        'realification_and_coordinate_order_checked': True,
        'primary_coefficients_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--primary-coefficients', type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    moments, matrices, support, contractions = generate()
    coefficient_data = {
        'schema': 'independent-semidefinite-lift-coefficients-v1',
        'coefficient_field': 'Q(sqrt(2))',
        'entry_encoding': '[row,column,rational_part,sqrt2_part], upper triangle, zero based',
        'dimension': 36,
        'moment_dimension': 47,
        'test_coordinate_groups': TEST_GROUPS,
        'test_coordinate_order': [f'{group}[{i}]' for group in TEST_GROUPS for i in range(4)],
        'moments': moments,
        'matrices': [
            {'index': k, 'moment': moments[k],
             'upper_entries': [[i, j, *value.json()] for (i, j), value in sorted(matrix.items())]}
            for k, matrix in enumerate(matrices)
        ],
    }
    coefficient_path = args.out_dir / 'independent_coefficients.json'
    json_dump(coefficient_path, coefficient_data)
    audit = redundancy_audit(moments, matrices)
    json_dump(args.out_dir / 'redundancy_certificate.json', audit)
    log = {
        'status': 'PASS',
        'arithmetic': 'Exact fractions in Q(sqrt(2),i), no floating point',
        'source_contractions_verified': contractions,
        'realification_order': 'Re first, Im first, Re second, Im second, u',
        'symmetric_lower_triangle_independently_expanded': True,
        'lower_triangle_complete_coefficient_dictionary_equality': True,
        'observed_laurent_support': support,
        'matrix_dimension': 36,
        'coefficient_matrix_count': 47,
        'nonzero_upper_entries_total': sum(len(matrix) for matrix in matrices),
        'coefficient_sha256': hashlib.sha256(coefficient_path.read_bytes()).hexdigest(),
        'common_kernel_dimension': audit['common_kernel_dimension'],
        'output_coefficient_family_common_vector_kernel_dimension': audit['output_coefficient_family_common_vector_kernel_dimension'],
        'auxiliary_coefficient_family_common_vector_kernel_dimension': audit['auxiliary_coefficient_family_common_vector_kernel_dimension'],
        'auxiliary_map_rank': audit['auxiliary_map_rank'],
        'full_coefficient_map_rank': audit['full_coefficient_map_rank'],
        'limitations': [
            'This program checks coefficients and the stated bounded redundancy audit.',
            'Exactness of the projected cone and all singular boundary claims require the mathematical proof.',
            'It does not check originality or global minimality.'
        ],
    }
    json_dump(args.out_dir / 'independent_verification_log.json', log)
    print(json.dumps(log, indent=2, sort_keys=True))
    if args.primary_coefficients is not None:
        comparison = compare_with_primary(args.primary_coefficients)
        json_dump(args.out_dir / 'independent_comparison.json', comparison)
        print(json.dumps(comparison, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
