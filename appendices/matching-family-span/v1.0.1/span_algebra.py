#!/usr/bin/env python3
"""Independent exact Laurent reconstruction of the pinned matching tests.

Python >= 3.10, standard library only.  This module never constructs the lift's
coefficient matrices and never performs the earlier redundancy audit.

Rhat = sqrt(2)*R0, so the returned test coordinates are
diag(sqrt(2)*I_32, I_4) times the authoritative real coefficient vector.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from math import lcm, prod


def rat(x):
    if not isinstance(x, (int, str, F)) or isinstance(x, bool):
        raise TypeError("Only exact integers, rational strings, and Fractions")
    return F(x)


@dataclass(frozen=True)
class G:
    """Gaussian rational, with no conversion through floating point."""
    re: F = F(0)
    im: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "re", rat(self.re))
        object.__setattr__(self, "im", rat(self.im))

    @staticmethod
    def of(x):
        return x if isinstance(x, G) else G(rat(x))

    def __add__(self, other):
        z = G.of(other)
        return G(self.re + z.re, self.im + z.im)
    __radd__ = __add__

    def __neg__(self):
        return G(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-G.of(other))

    def __rsub__(self, other):
        return G.of(other) - self

    def __mul__(self, other):
        z = G.of(other)
        return G(self.re*z.re - self.im*z.im,
                 self.re*z.im + self.im*z.re)
    __rmul__ = __mul__

    def __truediv__(self, other):
        z = G.of(other)
        d = z.norm()
        if not d:
            raise ZeroDivisionError
        p = self*z.conj()
        return G(p.re/d, p.im/d)

    def conj(self):
        return G(self.re, -self.im)

    def norm(self):
        return self.re*self.re + self.im*self.im

    def __bool__(self):
        return bool(self.re or self.im)


ZERO, ONE, I = G(), G(1), G(0, 1)


def clean(p):
    return {k: G.of(v) for k, v in p.items() if v}


def pc(x):
    return clean({0: G.of(x)})


def pa(p, q):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, ZERO) + v
    return clean(r)


def ps(p, z):
    return clean({k: v*z for k, v in p.items()})


def pm(p, q):
    r = {}
    for k, a in p.items():
        for ell, b in q.items():
            r[k+ell] = r.get(k+ell, ZERO) + a*b
    return clean(r)


def pj(p):
    """Conjugation on the unit circle: w^k -> w^(-k)."""
    return clean({-k: v.conj() for k, v in p.items()})


def pe(p, w):
    def power(k):
        z = ONE
        base = w if k >= 0 else w.conj()
        for _ in range(abs(k)):
            z = z*base
        return z
    return sum((v*power(k) for k, v in p.items()), ZERO)


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), ZERO)


def matching_vector(t):
    """Return (rational scaled coefficient vector, delta_hat, branch).

    This implements (12), (22), and (26) by formal Laurent matrix products;
    it does not use a pre-expanded four-coefficient formula or a rank search.
    """
    c, s, tau, rho = (rat(t[k]) for k in ("c", "s", "tau", "rho"))
    assert c*c + s*s == 1, "w0 is not on the unit circle"
    assert tau > 0, "The matching theorem requires tau > 0"
    u = [rat(x) for x in t["u"]]
    d = [[rat(x) for x in row] for row in t["D"]]
    assert len(u) == 4 and len(d) == 4 and all(len(row) == 4 for row in d)
    assert all(d[i][j] == d[j][i] for i in range(4) for j in range(4))
    w0 = G(c, s)
    wbar = w0.conj()
    # All entries are sqrt(2) times the original angular factor entries.
    r = [[-tau*wbar, -tau*w0],
         [tau*s*wbar*wbar + rho*I, tau*s*w0*w0 + rho*I],
         [tau*c*wbar*wbar + rho, -tau*c*w0*w0 + rho]]
    j = [[u[k+1]] + [-u[0] if k == ell else F(0) for ell in range(3)]
         for k in range(3)]
    sym = [[sum(j[a][k]*d[k][ell]*j[b][ell]
                for k in range(4) for ell in range(4))
            for b in range(3)] for a in range(3)]
    n = [[sum((sym[a][k]*r[k][b] for k in range(3)), ZERO)
          for b in range(2)] for a in range(3)]
    v = [n[2][b] - I*n[1][b] for b in range(2)]
    g, h = v[0] + v[1].conj(), v[0] - v[1].conj()
    delta = g.norm() + h.norm()
    if delta:
        assert delta > 0
        a = clean({0: G(g.norm()/delta), -2: h.norm()/delta*w0*w0})
        beta = clean({0: g*h/delta, -2: -g*h/delta*w0*w0})
        m = [[a, pj(beta)], [ps(beta, -1), pj(a)]]
        branch = "positive"
    else:
        assert not g and not h and all(not x for x in v)
        m = [[pc(1), {}], [{}, pc(1)]]
        branch = "zero"
    # Independently check the exact rotation, including the separate zero branch.
    for a in range(2):
        for b in range(2):
            entry = {}
            for k in range(2):
                entry = pa(entry, pm(pj(m[k][a]), m[k][b]))
            assert entry == pc(int(a == b)), "Rotation is not unitary"
            assert pe(m[a][b], w0) == G(int(a == b)), "M(alpha) != I"
    vm = [pa(ps(m[0][b], v[0]), ps(m[1][b], v[1])) for b in range(2)]
    assert pa(vm[0], pj(vm[1])) == pc(g)
    assert pa(vm[0], ps(pj(vm[1]), -1)) == clean({-2: h*w0*w0})
    # Chat = Rhat*M; Qtop_hat = -D*J^T*Chat.
    chat = [[pa(ps(m[0][b], r[a][0]), ps(m[1][b], r[a][1]))
             for b in range(2)] for a in range(3)]
    kreal = [[-sum(d[a][ell]*j[b][ell] for ell in range(4))
              for b in range(3)] for a in range(4)]
    top = [[{} for _ in range(2)] for _ in range(4)]
    for a in range(4):
        for b in range(2):
            for k in range(3):
                top[a][b] = pa(top[a][b], ps(chat[k][b], kreal[a][k]))
            assert set(top[a][b]) <= ({0, -2} if b == 0 else {0, 2})
    groups = [[top[a][b].get(power, ZERO) for a in range(4)]
              for b, power in ((0, 0), (0, -2), (1, 0), (1, 2))]
    vector = []
    for group in groups:
        vector.extend(z.re for z in group)
        vector.extend(z.im for z in group)
    vector.extend(u)
    assert len(vector) == 36
    return vector, delta, branch


def test_vector(t):
    if t["kind"] == "constant":
        index = t["index"]
        assert isinstance(index, int) and 0 <= index < 4
        vec = [F(0)]*36
        vec[index] = F(1)  # actual constant is e_index/sqrt(2)
        return vec, None, "constant"
    assert t["kind"] == "matching"
    return matching_vector(t)


def integerize_columns(matrix):
    """V -> W=V*diag(d_j) integral, with positive minimal d_j."""
    n = len(matrix)
    assert n and all(len(row) == n for row in matrix)
    scales = [lcm(*(matrix[i][j].denominator for i in range(n)))
              for j in range(n)]
    out = [[matrix[i][j]*scales[j] for j in range(n)] for i in range(n)]
    assert all(x.denominator == 1 for row in out for x in row)
    return [[int(x) for x in row] for row in out], scales


def bareiss(matrix):
    """Exact determinant certificate by integer fraction-free elimination."""
    a = [list(row) for row in matrix]
    n, previous, sign = len(a), 1, 1
    pivots, swaps = [], []
    for k in range(n-1):
        row = next((r for r in range(k, n) if a[r][k]), None)
        if row is None:
            return 0, {"pivots": pivots, "row_swaps": swaps}
        if row != k:
            a[k], a[row] = a[row], a[k]
            swaps.append([k, row])
            sign = -sign
        pivot = a[k][k]
        pivots.append(pivot)
        for i in range(k+1, n):
            for j in range(k+1, n):
                numerator = a[i][j]*pivot - a[i][k]*a[k][j]
                assert numerator % previous == 0, "Bareiss division is not exact"
                a[i][j] = numerator//previous
            a[i][k] = 0
        previous = pivot
    pivots.append(a[-1][-1])
    return sign*a[-1][-1], {"pivots": pivots, "row_swaps": swaps}


def determinant_certificate(matrix):
    integral, denominators = integerize_columns(matrix)
    det_int, trace = bareiss(integral)
    det_scaled = F(det_int, prod(denominators))
    return {
        "column_denominators": denominators,
        "integer_matrix": integral,
        "integer_determinant": det_int,
        "scaled_determinant": str(det_scaled),
        "actual_determinant": str(det_scaled/F(2**16)),
        "bareiss": trace,
    }
