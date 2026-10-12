#!/usr/bin/env python3
"""Independent exact checks for the matching and singular-boundary audit.

Dependency: SymPy 1.14.0 (other recent versions should also work).
Run from any directory: python check_matching.py

All identities are polynomial/Laurent/rational identities, with no floating
point arithmetic.  This checker intentionally does not import the lift's
coefficient generator.  Generic real symmetric symbols never impose positive
semidefiniteness.  The analytic implication from the scalar certificates to
cone membership is proved in matching_boundary_audit.md.
"""


if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

import json
import sympy as s


I = s.I
rt2 = s.sqrt(2)
w, w0 = s.symbols("w w0", nonzero=True)
conjugates = {s.conjugate(w): 1 / w, s.conjugate(w0): 1 / w0}
results = []


def complex_symbols(*names):
    result = []
    for name in names:
        a, b = s.symbols(name + " " + name + "_bar")
        conjugates[s.conjugate(a)] = b
        conjugates[s.conjugate(b)] = a
        result.append(a)
    return result


def star(a):
    return s.conjugate(a).xreplace(conjugates)


def adj(a):
    return a.applyfunc(star).T


def real(a):
    return (a + star(a)) / 2


def normal(a):
    return s.cancel(s.expand(a))


def zero(a):
    entries = list(a) if isinstance(a, s.MatrixBase) else [a]
    for entry in entries:
        assert not entry.atoms(s.Float), "Floating point value encountered"
        assert normal(entry) == 0, str(normal(entry))


def check(name, a):
    zero(a)
    results.append(name)


def frequencies(a):
    numerator, denominator = normal(a).as_numer_denom()
    if numerator == 0:
        return set()
    denominator_terms = s.Poly(denominator, w).terms()
    assert len(denominator_terms) == 1, "Expression is not Laurent in w"
    denominator_power = denominator_terms[0][0][0]
    return {powers[0]-denominator_power
            for powers, coefficient in s.Poly(numerator, w).terms()
            if normal(coefficient) != 0}


def U(t):
    c = (t + 1 / t) / 2
    sn = (t - 1 / t) / (2 * I)
    return s.Matrix([[-1 / t, -t], [sn / t**2, sn * t**2],
                     [c / t**2, -c * t**2]]) / rt2


V = s.Matrix([[0, 0], [I, I], [1, 1]]) / rt2


def Q(y):
    a, b, c = y
    return s.Matrix([[a*a + b*b, -a*b, -a*c],
                     [-a*b, b*b + c*c, -b*c],
                     [-a*c, -b*c, c*c + a*a]])


def polar(y, z):
    return (Q(y + z) - Q(y) - Q(z)) / 2


def pairing(a, b, S):
    return real(s.trace(adj(a) * S * b))


def B(r):
    return s.Matrix([[0, r[0], r[1], r[2]],
                     [-r[0], 0, 0, 0], [-r[1], 0, 0, 0],
                     [-r[2], 0, 0, 0]])


def J(u):
    return s.Matrix.hstack(s.Matrix(u[1:]), -u[0] * s.eye(3))


def Phi_from_Q(W, q):
    return s.Matrix.vstack(
        s.Matrix.hstack(s.Matrix([[s.trace(q * W[1:, 1:])]]),
                        -W[0:1, 1:] * q),
        s.Matrix.hstack(-q * W[1:, 0:1], W[0, 0] * q))


def realify(F):
    columns = []
    for j in range(2):
        columns.extend([F[:, j].applyfunc(real),
                        F[:, j].applyfunc(lambda x: (x-star(x))/(2*I))])
    return s.Matrix.hstack(*columns)


# 1. Universal angular Gram identity, using the exact source order.
eta, xi = s.symbols("eta xi", real=True)
c, sn = (w + 1 / w) / 2, (w - 1 / w) / (2 * I)
y = s.Matrix([eta*c, eta*sn, xi])
F = eta * U(w) + xi * V
RF = realify(F)
check("realification_order_and_angular_Gram", RF * RF.T - Q(y))
check("realification_equals_real_complex_Gram",
      RF * RF.T - (F * adj(F)).applyfunc(real))

# 2. Universal B, J, Phi identities: one rank-one Gram summand suffices,
# because all identities are linear in the 3x3 Gram matrix.
r = s.Matrix(s.symbols("r1:4", real=True))
u = s.Matrix(s.symbols("u0:4", real=True))
W = s.zeros(4)
for row in range(4):
    for col in range(row, 4):
        W[row, col] = W[col, row] = s.Symbol(f"W{row}{col}", real=True)
check("B_r_u_equals_J_u_transpose_r", B(r)*u - J(u).T*r)
check("Phi_rank_one_Gram_summand", B(r).T*W*B(r) - Phi_from_Q(W, r*r.T))
check("Phi_scalar_trace_identity",
      (u.T*Phi_from_Q(W, r*r.T)*u)[0] - s.trace(J(u)*W*J(u).T*r*r.T))
beta = s.Matrix(s.symbols("beta1:4", real=True))
Xgeneric = s.zeros(4)
for row in range(4):
    for col in range(row, 4):
        Xgeneric[row, col] = Xgeneric[col, row] = s.Symbol(f"X{row}{col}", real=True)
qtop = -W*J(u).T*beta  # Here generic symmetric W plays the role of D.
check("generic_indefinite_D_certificate_top_summand",
      (qtop.T*Xgeneric*qtop)[0] - s.trace(J(u)*W*Xgeneric*W*J(u).T*beta*beta.T))
check("generic_indefinite_D_certificate_cross_summand",
      2*(qtop.T*B(r)*u)[0] + 2*(beta.T*J(u)*W*J(u).T*r)[0])

# 3. Universal rotation identities, with g,h unconstrained complex numbers.
g, h = complex_symbols("g", "h")
gb, hb = star(g), star(h)
delta = g*gb + h*hb
kappa = (w0/w)**2
a = (g*gb + h*hb*kappa) / delta
bbar = g*h*(1-kappa) / delta
b = star(bbar)
M = s.Matrix([[a, b], [-bbar, star(a)]])
A = s.Matrix([[g, hb], [h, -gb]])
check("rotation_coefficient_matrix_scaled_unitary", adj(A)*A - delta*s.eye(2))
check("rotation_linear_system", A*s.Matrix([a, bbar]) - s.Matrix([g, h*kappa]))
check("rotation_unitarity", M*adj(M) - s.eye(2))
check("rotation_at_base_angle", M.subs(w, w0) - s.eye(2))
v = s.Matrix([[(g+h)/2, (gb-hb)/2]])
vM = v*M
check("rotation_sum_pairing", vM[0] + star(vM[1]) - g)
check("rotation_difference_pairing", vM[0] - star(vM[1]) - h*kappa)
rp, rm = complex_symbols("rp", "rm")
rowM = s.Matrix([[rp, rm]])*M
assert frequencies(rowM[0]) <= {0, -2}
assert frequencies(rowM[1]) <= {0, 2}
results.append("every_C_first_column_0_minus2_second_0_plus2")
assert frequencies(star(rowM[1])) <= {0, -2}
results.append("conjugated_second_column_frequency_bound")

# 4. Independently verify the displayed cubic-frequency expansion for an
# arbitrary complex 3x2 matrix, before using the rotation constraints.
T = s.Matrix(3, 2, complex_symbols("T11", "T12", "T21", "T22", "T31", "T32"))
tv = T[2, :] - I*T[1, :]
tm = T[2, :] + I*T[1, :]
tn = T[0, :]
reduced = real(w**3*(tv[0] - star(tv[1]))/2 +
               w*(tm[0]/2-tn[0]-star(tm[1])/2-star(tn[1])))
check("universal_U_pairing_frequency_expansion",
      rt2*pairing(U(w), T, s.eye(3)) - reduced)
check("universal_V_pairing_expansion",
      rt2*pairing(V, T, s.eye(3)) - real(tv[0]+tv[1]))

# 5. Generic symmetric S: no sign conditions or definiteness assumptions.
S = s.zeros(3)
for row in range(3):
    for col in range(row, 3):
        S[row, col] = S[col, row] = s.Symbol(f"S{row}{col}", real=True)
tau, rho = s.symbols("tau rho", real=True)
ca, sa = (w0+1/w0)/2, (w0-1/w0)/(2*I)
da, ka = s.Matrix([ca, sa, 0]), s.Matrix([0, 0, 1])
z = tau*da + rho*ka
R0 = tau*U(w0) + rho*V
check("generic_S_V_pairing_at_base_angle", pairing(V, R0, S) - s.trace(S*polar(ka, z)))
check("generic_S_U_pairing_at_base_angle", pairing(U(w0), R0, S) - s.trace(S*polar(da, z)))
Up = U(w0).diff(w0) * I*w0
check("generic_S_U_derivative_at_base_angle",
      pairing(Up, R0, S) - s.trace(S*polar(s.Matrix([-sa, ca, 0]), z)))
Mprime = (M.diff(w)*I*w).subs(w, w0)
check("rotation_derivative_is_skew_Hermitian", Mprime + adj(Mprime))
g0, g1, gr, gi = s.symbols("G0 G1 Gr Gi", real=True)
G = s.Matrix([[g0, gr+I*gi], [gr-I*gi, g1]])
check("indefinite_Hermitian_stationarity", real(s.trace(G*Mprime)))

# 6. Delta=0 is checked as its own algebraic branch, without division by 0.
# g=h=0 implies v=0; hence N3=i*N2, and C=R0 (M=I).
n1, n2, n3, n4 = complex_symbols("n1", "n2", "n3", "n4")
Nzero = s.Matrix([[n1, n2], [n3, n4], [I*n3, I*n4]])
check("delta_zero_V_pairing", pairing(V, Nzero, s.eye(3)))
assert frequencies(pairing(U(w), Nzero, s.eye(3))) <= {-1, 1}
results.append("delta_zero_U_pairing_only_first_harmonic")
Szero = s.diag(1, -1, 0)  # Nonzero and indefinite, with z=(1,0,0).
Rzero = U(s.Integer(1))
check("delta_zero_nonzero_indefinite_S_example_v",
      (Szero*Rzero)[2, :] - I*(Szero*Rzero)[1, :])
check("delta_zero_nonzero_indefinite_S_full_matching",
      pairing(F, Rzero, Szero) - s.trace(Szero*polar(y, s.Matrix([1, 0, 0]))))

# 7. One exact, nonzero-delta, indefinite-D combined matching regression.
# The universal proof is the preceding symbolic decomposition, not this test.
ue = s.Matrix([1, 2, -1, 3])
De = s.diag(-2, 1, 3, -1)
Se = J(ue)*De*J(ue).T
assert Se[0, 0] < 0 and Se[1, 1] > 0
results.append("combined_example_D_and_S_are_indefinite")
Re = U(s.Integer(1)) + 2*V
Ne = Se*Re
ve = Ne[2, :] - I*Ne[1, :]
ge, he = normal(ve[0]+star(ve[1])), normal(ve[0]-star(ve[1]))
de = normal(ge*star(ge)+he*star(he))
assert de > 0
Me = M.subs({g:ge, gb:star(ge), h:he, hb:star(he), w0:1}, simultaneous=True)
Ce = Re*Me
check("indefinite_D_combined_matching_Gram", (Ce*adj(Ce)).applyfunc(real)-Q(s.Matrix([1,0,2])))
check("indefinite_D_combined_matching_pairing",
      pairing(F, Ce, Se)-s.trace(Se*polar(y, s.Matrix([1,0,2]))))

# 8. Exact singular-boundary witnesses, including the z1=z2=0 axis.
# These verify the completion and every B_j range relation in examples;
# the quantified boundary argument is analytic in the companion audit note.
for label, ys, Xs, Fs in [
        ("rank3_nonaxis", s.Matrix([1,0,0]), s.diag(1,1,0,1), U(s.Integer(1))),
        ("rank3_axis", s.Matrix([0,0,1]), s.diag(1,0,1,1), V)]:
    Rs = realify(Fs)
    Bs = [B(Rs[:, j]) for j in range(4)]
    for j, Bj in enumerate(Bs, start=1):
        check(f"{label}_range_B{j}", Xs*Bj-Bj)
    Zs = Phi_from_Q(Xs, Q(ys))
    stacked = s.Matrix.vstack(*Bs)
    top = s.kronecker_product(s.eye(4), Xs)
    Ls = s.Matrix.vstack(s.Matrix.hstack(top, stacked), s.Matrix.hstack(stacked.T, Zs))
    factor = s.Matrix.hstack(top, stacked)
    check(f"{label}_exact_PSD_Gram_completion", Ls-factor.T*factor)
    assert Xs.rank() == 3 and Ls.rank() == 12
    results.append(f"{label}_ranks_X3_L12")

# 9. The constant I4/I4 output coefficient combination is I36, in the
# independent note's explicit real basis. This is an exact common-kernel
# certificate, rather than a floating point rank computation.
def basis_at(c2, s2):
    P = s.zeros(20, 36)
    for k in range(4):
        P[k, k] = 1
        P[4+k, 4+k] = 1
        P[k, 8+k], P[4+k, 8+k] = c2, -s2
        P[k, 12+k], P[4+k, 12+k] = s2, c2
        P[8+k, 16+k] = 1
        P[12+k, 20+k] = 1
        P[8+k, 24+k], P[12+k, 24+k] = c2, s2
        P[8+k, 28+k], P[12+k, 28+k] = -s2, c2
        P[16+k, 32+k] = 1
    return P


P0, P1 = basis_at(1, 0), basis_at(-1, 0)
check("common_kernel_identity_combination_I36", (P0.T*P0 + P1.T*P1)/2-s.eye(36))

print(json.dumps({"status": "PASS", "arithmetic": "exact symbolic; no floats",
                  "sympy_version": s.__version__, "checks": results,
                  "check_count": len(results),
                  "limitations": ["Not a formal proof assistant certificate",
                                  "Cone equivalence uses the accompanying analytic proof",
                                  "The sample regressions supplement universal identities"]}, indent=2))
