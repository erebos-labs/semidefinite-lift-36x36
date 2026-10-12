#!/usr/bin/env python3
"""Exact denominator-free block-elimination certificate for det(L)=p.

Requires SymPy. The determinant argument following the checked matrix
identity is recorded explicitly in the JSON output and THEOREM_PROOF.md.
It cancels a nonzero polynomial in a polynomial ring, never assumes that
the evaluated X is invertible, and never changes the degree-20 p.
"""

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

import json
import sympy as s


def symmetric(n, prefix):
    m=s.zeros(n)
    for i in range(n):
        for j in range(i,n):
            m[i,j]=m[j,i]=s.Symbol(f"{prefix}{i}{j}")
    return m


def check_zero(m):
    for value in m:
        assert not value.atoms(s.Float)
        assert s.expand(value)==0


def B(r):
    return s.Matrix([[0,r[0],r[1],r[2]],[-r[0],0,0,0],
                     [-r[1],0,0,0],[-r[2],0,0,0]])


x,z=symmetric(4,"x"),symmetric(4,"z")
d=x.det(method="domain-ge")
adj=x.adjugate(method="domain-ge")
check_zero(x*adj-d*s.eye(4))
check_zero(adj*x-d*s.eye(4))
assert d.subs({x[i,j]:int(i==j) for i in range(4) for j in range(i,4)})==1
r=[s.Matrix([s.Symbol(f"r{j}_{k}") for k in range(3)]) for j in range(4)]
b=s.Matrix.vstack(*[B(v) for v in r])
gram=sum((v*v.T for v in r),s.zeros(3))
phi=s.Matrix.vstack(
    s.Matrix.hstack(s.Matrix([[s.trace(gram*adj[1:,1:])]]),-adj[0:1,1:]*gram),
    s.Matrix.hstack(-gram*adj[1:,0:1],adj[0,0]*gram))
top=s.kronecker_product(s.eye(4),x)
hat_adj=s.kronecker_product(s.eye(4),adj)
check_zero(b.T*hat_adj*b-phi)
l=s.Matrix.vstack(s.Matrix.hstack(top,b),s.Matrix.hstack(b.T,z))
right=s.Matrix.vstack(s.Matrix.hstack(d*s.eye(16),-hat_adj*b),
                      s.Matrix.hstack(s.zeros(4,16),d*s.eye(4)))
target=s.Matrix.vstack(s.Matrix.hstack(d*top,s.zeros(16,4)),
                       s.Matrix.hstack(d*b.T,d*z-phi))
check_zero(l*right-target)
print(json.dumps({
    "status":"PASS", "arithmetic":"exact polynomial, no division or floats",
    "sympy_version":s.__version__,
    "checks":["X adj(X)=det(X) I", "adj(X) X=det(X) I",
              "det(X) is a nonzero polynomial (value 1 at I)",
              "four arbitrary skew blocks give Phi of their Gram",
              "all 400 entries of L R equal the displayed block-triangular target"],
    "determinant_deduction":{
        "det_R":"d^20",
        "det_upper_left_target":"det(d (I4 tensor X)) = d^16 d^4 = d^20",
        "det_lower_right_target":"det(d Z-Phi_y(adj X)) = the original p",
        "identity":"d^20 det(L)=d^20 p",
        "cancellation":"Cancel the nonzero polynomial d^20 in the integral polynomial ring",
        "conclusion":"det(L)=p identically, hence for singular and nonsingular X",
        "degree_reduction":"none; p remains homogeneous of degree 20"},
    "limits":["Uses the angular Gram identity checked separately",
              "Does not expand the full degree-20 polynomial monomial by monomial",
              "Determinant multiplicativity, triangular determinants, and integral-domain cancellation are mathematical proof rules, not Lean-checked"]
},indent=2))
