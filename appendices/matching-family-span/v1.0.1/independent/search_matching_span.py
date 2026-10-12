#!/usr/bin/env python3
"""One predeclared finite exact search for a rank-36 matching-test witness.

Source: erebos-labs/semidefinite-lift-36x36 at commit
2d2c170206ab937bcf56ff3290066cc3b2dfdd99, THEOREM_PROOF.md §§3,6,8.

This does not reconstruct the old lift coefficients, repeat a redundancy
audit, or search alternative lifts.  It constructs only the admissible
matching tests in equation (26), using exactly the two branches of §6,
and four constant tests in the first four-vector block.

The candidate grid is fixed below before any rank computation.  It has
1,881 candidates, tested in w0,rho,u,D order. Search stops immediately upon
rank 36. No adaptive grid expansion is implemented.

All arithmetic is Fraction, using Gaussian rational pairs.  Write
Rhat=sqrt(2)*R0.  The formula for M is unchanged when g,h are both scaled
by sqrt(2). The first 32 coordinates of the resulting q have a common
1/sqrt(2) factor. Scaling those rows by sqrt(2) makes each witness column
rational. The constant tests are e_i/sqrt(2), i=0,...,3, so they obey the
same row scaling.  Rank is computed by exact rational elimination.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path


@dataclass(frozen=True)
class GQ:
    """A Gaussian rational number, with no floating point conversion."""
    re: F = F(0)
    im: F = F(0)

    def __add__(self, other):
        if not isinstance(other, GQ):
            other = GQ(F(other))
        return GQ(self.re+other.re, self.im+other.im)

    __radd__ = __add__

    def __neg__(self):
        return GQ(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-other if isinstance(other, GQ) else -GQ(F(other)))

    def __mul__(self, other):
        if not isinstance(other, GQ):
            other = GQ(F(other))
        return GQ(self.re*other.re-self.im*other.im,
                  self.re*other.im+self.im*other.re)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, GQ):
            other = GQ(F(other))
        denominator = other.re**2+other.im**2
        if denominator == 0:
            raise ZeroDivisionError
        numerator = self*other.conjugate()
        return GQ(numerator.re/denominator, numerator.im/denominator)

    def conjugate(self):
        return GQ(self.re, -self.im)

    def norm2(self):
        return self.re**2+self.im**2

    def json(self):
        return [str(self.re), str(self.im)]


ZERO = GQ()
I = GQ(F(0), F(1))
PIN = '2d2c170206ab937bcf56ff3290066cc3b2dfdd99'
TEST_LABELS = [f'{letter}.{component}.{k}'
               for letter in ['A', 'B', 'C', 'D']
               for component in ['real', 'imag'] for k in range(4)] + [f'u.{k}' for k in range(4)]


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(column) for column in zip(*a)]


def matmul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), 0)
             for j in range(len(b[0]))] for i in range(len(a))]


def fixed_grid():
    """The complete declared parameter grid, independent of rank results."""
    angles = [GQ(F(1)), GQ(F(0), F(1)), GQ(F(3, 5), F(4, 5))]
    rhos = [F(0), F(1), F(-1)]
    us = []
    for i in range(4):
        us.append((f'e{i}', [F(k == i) for k in range(4)]))
    for i in range(4):
        for j in range(i+1, 4):
            us.append((f'e{i}+e{j}', [F(k in [i, j]) for k in range(4)]))
    us.append(('ones', [F(1)]*4))
    matrices = [('zero', [[F(0)]*4 for _ in range(4)]), ('identity', identity(4))]
    for i in range(4):
        matrices.append((f'E{i}{i}', [[F(a == i and b == i) for b in range(4)] for a in range(4)]))
    for i in range(4):
        for j in range(i+1, 4):
            matrices.append((f'F{i}{j}', [[F((a, b) in [(i, j), (j, i)]) for b in range(4)] for a in range(4)]))
    for i in range(4):
        for j in range(i+1, 4):
            matrices.append((f'I+F{i}{j}', [[F(a == b)+F((a, b) in [(i, j), (j, i)]) for b in range(4)] for a in range(4)]))
    matrices.append(('alternating', [[F([1, -1, 1, -1][a]) if a == b else F(0) for b in range(4)] for a in range(4)]))
    assert len(angles) == 3 and len(rhos) == 3 and len(us) == 11 and len(matrices) == 19
    return angles, rhos, us, matrices


def grid_description():
    angles, rhos, us, matrices = fixed_grid()
    return {
        'source_repository': 'https://github.com/erebos-labs/semidefinite-lift-36x36',
        'source_commit': PIN,
        'tau': '1',
        'angles': [{'c': str(w.re), 's': str(w.im)} for w in angles],
        'rho_values': [str(rho) for rho in rhos],
        'u_values': [{'label': name, 'u': list(map(str, u))} for name, u in us],
        'D_values': [{'label': name, 'D': [list(map(str, row)) for row in D]} for name, D in matrices],
        'iteration_order': ['angle', 'rho', 'u', 'D'],
        'maximum_matching_candidates': len(angles)*len(rhos)*len(us)*len(matrices),
        'initial_constants': [{'kind': 'constant', 'index': i, 'normalization': '1/sqrt(2)'} for i in range(4)],
        'stopping_rule': 'Stop immediately at exact rank36; otherwise stop after all1881 candidates. No grid expansion.',
    }


def matching_column(c, s, tau, rho, u, D):
    """Compact exact coefficients of the authoritative formulas (22),(26).

    Returns the rational column after scaling the first32 rows by sqrt(2),
    plus branch metadata.  Rhat and ghat,hhat are scaled by sqrt(2), so
    delta_hat=2*delta.  The branch condition is unaffected.
    """
    c, s, tau, rho = map(F, (c, s, tau, rho))
    u = list(map(F, u))
    D = [list(map(F, row)) for row in D]
    assert tau > 0 and c*c+s*s == 1
    assert len(u) == 4 and len(D) == 4 and all(len(row) == 4 for row in D)
    assert all(D[i][j] == D[j][i] for i in range(4) for j in range(4))
    w = GQ(c, s)
    wc = w.conjugate()
    w2, wc2 = w*w, wc*wc
    Uhat = [[-wc, -w], [s*wc2, s*w2], [c*wc2, -c*w2]]
    Vhat = [[ZERO, ZERO], [I, I], [GQ(F(1)), GQ(F(1))]]
    Rhat = [[tau*Uhat[i][j]+rho*Vhat[i][j] for j in range(2)] for i in range(3)]
    J = [[u[i+1]]+[F(-u[0] if i == j else 0) for j in range(3)] for i in range(3)]
    DJT = matmul(D, transpose(J))
    S = matmul(J, DJT)
    Nhat = matmul(S, Rhat)
    v = [Nhat[2][j]-I*Nhat[1][j] for j in range(2)]
    g, h = v[0]+v[1].conjugate(), v[0]-v[1].conjugate()
    delta_hat = g.norm2()+h.norm2()
    assert delta_hat >= 0
    if delta_hat:
        branch = 'delta>0'
        a_constant = GQ(g.norm2()/delta_hat)
        a_minus2 = (h.norm2()/delta_hat)*w2
        bbar_constant = g*h/delta_hat
        bbar_minus2 = -bbar_constant*w2
        C1_constant = [Rhat[i][0]*a_constant-Rhat[i][1]*bbar_constant for i in range(3)]
        C1_minus2 = [Rhat[i][0]*a_minus2-Rhat[i][1]*bbar_minus2 for i in range(3)]
        C2_constant = [Rhat[i][0]*bbar_constant.conjugate()+Rhat[i][1]*a_constant.conjugate() for i in range(3)]
        C2_plus2 = [Rhat[i][0]*bbar_minus2.conjugate()+Rhat[i][1]*a_minus2.conjugate() for i in range(3)]
    else:
        branch = 'delta=0'
        assert not g.norm2() and not h.norm2()
        # The source requires the direct branch M=I, not a limit.
        C1_constant = [Rhat[i][0] for i in range(3)]
        C1_minus2 = [ZERO]*3
        C2_constant = [Rhat[i][1] for i in range(3)]
        C2_plus2 = [ZERO]*3
    K = [[-x for x in row] for row in DJT]
    groups = []
    for column in [C1_constant, C1_minus2, C2_constant, C2_plus2]:
        groups.append([sum((K[i][j]*column[j] for j in range(3)), ZERO) for i in range(4)])
    scaled_column = []
    for vector in groups:
        scaled_column.extend(z.re for z in vector)
        scaled_column.extend(z.im for z in vector)
    scaled_column.extend(u)
    assert len(scaled_column) == 36
    return scaled_column, {
        'delta_hat': str(delta_hat), 'delta': str(delta_hat/2), 'branch': branch,
        'g_hat': g.json(), 'h_hat': h.json(),
    }


class ExactSpan:
    def __init__(self):
        self.pivots = {}

    def add(self, column):
        vector = list(column)
        for pivot in sorted(self.pivots):
            value = vector[pivot]
            if value:
                row = self.pivots[pivot]
                vector = [x-value*y for x, y in zip(vector, row)]
        pivot = next((k for k, x in enumerate(vector) if x), None)
        if pivot is None:
            return None
        value = vector[pivot]
        self.pivots[pivot] = [x/value for x in vector]
        return pivot, value

    def rank(self):
        return len(self.pivots)


def determinant(rows):
    a = [row[:] for row in rows]
    n, result = len(a), F(1)
    assert all(len(row) == n for row in a)
    for column in range(n):
        pivot = next((r for r in range(column, n) if a[r][column]), None)
        if pivot is None:
            return F(0)
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            result = -result
        value = a[column][column]
        result *= value
        for row in range(column+1, n):
            if a[row][column]:
                multiple = a[row][column]/value
                for k in range(column, n):
                    a[row][k] -= multiple*a[column][k]
    return result


def search():
    span = ExactSpan()
    selected, columns = [], []
    for i in range(4):
        column = [F(k == i) for k in range(36)]
        pivot, value = span.add(column)
        selected.append({'kind': 'constant', 'index': i,
                         'normalization': '1/sqrt(2)',
                         'selected_column': len(selected),
                         'rank_after': span.rank(), 'pivot_coordinate': pivot,
                         'pivot_value': str(value)})
        columns.append(column)
    angles, rhos, us, matrices = fixed_grid()
    scanned = 0
    branches_scanned = {'delta=0': 0, 'delta>0': 0}
    for angle_index, w in enumerate(angles):
        for rho_index, rho in enumerate(rhos):
            for u_index, (u_label, u) in enumerate(us):
                for D_index, (D_label, D) in enumerate(matrices):
                    candidate_index = scanned
                    scanned += 1
                    column, metadata = matching_column(w.re, w.im, F(1), rho, u, D)
                    branches_scanned[metadata['branch']] += 1
                    new_pivot = span.add(column)
                    if new_pivot is not None:
                        pivot, value = new_pivot
                        selected.append({
                            'kind': 'matching', 'candidate_index': candidate_index,
                            'grid_indices': [angle_index, rho_index, u_index, D_index],
                            'c': str(w.re), 's': str(w.im), 'tau': '1', 'rho': str(rho),
                            'u': list(map(str, u)), 'D': [list(map(str, row)) for row in D],
                            'u_label': u_label, 'D_label': D_label,
                            'selected_column': len(selected),
                            'rank_after': span.rank(), 'pivot_coordinate': pivot,
                            'pivot_value': str(value), **metadata,
                        })
                        columns.append(column)
                    if span.rank() == 36:
                        return span, selected, columns, scanned, branches_scanned
    return span, selected, columns, scanned, branches_scanned


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def main():
    if not __debug__:
        raise RuntimeError('Run without -O: exact assertions must remain enabled.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    grid_path = args.out_dir/'declared_grid.json'
    # Preserve the declared finite grid before running the search.
    write_json(grid_path, grid_description())
    span, selected, columns, scanned, branches_scanned = search()
    rows = transpose(columns)
    det = determinant(rows) if span.rank() == 36 else None
    if det is not None:
        assert det != 0
    branches_selected = {branch: sum(test.get('branch') == branch for test in selected)
                         for branch in ['delta=0', 'delta>0']}
    witness = {
        'schema': 'matching-test-span-witness-v1',
        'source_commit': PIN,
        'source_equations': ['THEOREM_PROOF.md section3 equation6',
                             'THEOREM_PROOF.md section6 equation22',
                             'THEOREM_PROOF.md section8 equation26'],
        'tests': selected,
        'rank': span.rank(),
        'test_coordinate_labels': TEST_LABELS,
        'scaled_matrix_rows': [list(map(str, row)) for row in rows],
        'scaled_columns': [list(map(str, col)) for col in columns],
        'row_scaling': 'scaled=diag(sqrt(2)*I32,I4)*original; constants are e_i/sqrt(2)',
        'original_entry_rule': 'For row<32, value=scaled_entry*sqrt(2)/2; otherwise value=scaled_entry.',
        'scaled_matrix_determinant': str(det) if det is not None else None,
        'original_matrix_determinant': str(det/F(2**16)) if det is not None else None,
        'branches_selected': branches_selected,
    }
    witness_path = args.out_dir/'span_witness.json'
    write_json(witness_path, witness)
    log = {
        'status': 'FULL_RANK_WITNESS' if span.rank() == 36 else 'GRID_EXHAUSTED_NO_FULL_RANK_WITNESS',
        'rank': span.rank(), 'constants_selected': 4,
        'matching_tests_selected': len(selected)-4,
        'matching_candidates_scanned': scanned,
        'maximum_matching_candidates': 1881,
        'stopped_immediately_at_full_rank': span.rank() == 36,
        'branches_scanned': branches_scanned,
        'branches_selected': branches_selected,
        'scaled_matrix_determinant': witness['scaled_matrix_determinant'],
        'original_matrix_determinant': witness['original_matrix_determinant'],
        'declared_grid_sha256': hashlib.sha256(grid_path.read_bytes()).hexdigest(),
        'witness_sha256': hashlib.sha256(witness_path.read_bytes()).hexdigest(),
        'limitations': [
            'One finite grid; no adaptive expansion and no alternative lift search.',
            'This determines the span of the exhibited admissible tests, not lower bounds on arbitrary semidefinite lifts.',
            'Source matching admissibility is reconstructed from exact formulas; a separate verifier should check this witness independently.'
        ],
    }
    write_json(args.out_dir/'search_log.json', log)
    print(json.dumps(log, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
