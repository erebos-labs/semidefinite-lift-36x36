# Independent matching and boundary audit

## Scope and outcome

This audit independently read `source/build/main.tex`, all five requested
section files, and `source/build/references.bib` from the supplied snapshot of
OpenAI's *An Exact Semidefinite Lift of a Nonspectrahedral Hyperbolicity Cone*,
at commit `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`.

**Outcome:** the matching argument in the pinned source does fit the proposed
36-dimensional test space. Its converse proof continues to work on the entire
closed degree-20 cone with the restricted 47-dimensional functional domain.
No extra assumption on the sign of the matching matrix, no representing
measure, no continuity of matching factors, and no positive-definite
approximation of a singular matrix is needed. This audit found no failure of
the proposed 36/24 construction. The primary deliverable's coefficient
generator and auxiliary-rank computation are separate verification components;
this note does not substitute for them.

The source itself states a homogeneous real symmetric size-100 pencil with
307 auxiliary variables. Its calculation is explicit: a 20-dimensional slice
pencil tested against five scalar harmonics gives size 100; 22 independent
slice coordinates times 15 harmonics give 330 functional coordinates, of
which 23 are prescribed outputs. This verifies the claimed published baseline
directly from `introduction.tex` and `lift.tex`.

The source's matching identity is inherited mathematics. Detecting its
one-sided column supports and pruning unused moment coordinates are the
compression used here. Neither this note nor its checker establishes
historical originality or a global lower bound on lift size.

## 1. Fixed polynomial and the exact closed slice

Use the source's exact matrix

\[
Q(y)=\begin{pmatrix}
y_1^2+y_2^2&-y_1y_2&-y_1y_3\\
-y_1y_2&y_2^2+y_3^2&-y_2y_3\\
-y_1y_3&-y_2y_3&y_3^2+y_1^2
\end{pmatrix}
\]

and

\[
\Phi_y\begin{pmatrix}a&r^T\\r&T\end{pmatrix}
=\begin{pmatrix}\operatorname{tr}(Q(y)T)&-r^TQ(y)\\
-Q(y)r&aQ(y)\end{pmatrix}.
\]

The defining polynomial is
\(p(X,Z,y)=\det((\det X)Z-\Phi_y(\operatorname{adj}X))\), of homogeneous
degree 20, and \(e=(I_4,I_4,0)\). Cone membership means that **all 20 roots,
with multiplicity**, of \(p(te-(X,Z,y))\) are nonnegative. Zero roots are
retained. No degree-16 polynomial is substituted.

Put \(w=e^{i\theta}\), \(c=\cos\theta\), \(s=\sin\theta\), and

\[
U(\theta)=\frac1{\sqrt2}\begin{pmatrix}
-w^{-1}&-w\\sw^{-2}&sw^2\\cw^{-2}&-cw^2
\end{pmatrix},\qquad
V=\frac1{\sqrt2}\begin{pmatrix}0&0\\i&i\\1&1\end{pmatrix}.
\]

For \(F=\eta U+\xi V=(f_1,f_2)\), its realification is exactly
\(\widetilde F=(\Re f_1,\Im f_1,\Re f_2,\Im f_2)\).
Exact multiplication gives
\(\widetilde F\widetilde F^T=Q(\eta c,\eta s,\xi)\). The independent
checker verifies this universal identity as a Laurent-polynomial identity.

For \(r\in\mathbb R^3\), put

\[
B(r)=\begin{pmatrix}0&r^T\\-r&0\end{pmatrix},\qquad
J_u=(u'\; -u_0I_3),\quad u=(u_0,u').
\]

The checker verifies \(B(r)u=J_u^Tr\), and verifies, with a completely
generic real symmetric \(W\),
\(B(r)^TWB(r)=\Phi_{rr^T}(W)\), where the latter notation means replacing
\(Q(y)\) by \(rr^T\) in the displayed block formula. Adding the four
rank-one identities gives, for the ordered columns \(r_j\) of
\(\widetilde F\),

\[
\sum_{j=1}^4B(r_j)^TWB(r_j)=\Phi_y(W),\qquad
u^T\Phi_y(W)u=\operatorname{tr}(J_uWJ_u^TQ(y)).
\]

Put \(\mathcal B=\operatorname{col}_{j=1}^4 B(r_j)\). The size-20 pencil is

\[
L_\theta=\begin{pmatrix}
I_4\otimes\mathcal X&\mathcal B\\
\mathcal B^T&\mathcal Z
\end{pmatrix},
\]

For invertible \(\mathcal X\), the block determinant identity gives

\[
\det L_\theta=(\det\mathcal X)^4
\det(\mathcal Z-\Phi_y(\mathcal X^{-1}))
=p(\mathcal X,\mathcal Z,y).
\]

Both ends are polynomials. Their identity on the nonempty open set of
invertible \(\mathcal X\) is a polynomial identity, hence holds at every
singular \(\mathcal X\) too. Since \(L_\theta(I,I,0,0)=I_{20}\), linearity
then gives the exact identity

\[
p(te-(\mathcal X,\mathcal Z,Y_\theta(\eta,\xi)))
=\det(tI_{20}-L_\theta(\mathcal X,\mathcal Z,\eta,\xi)).
\]

Every \(y\) lies on such a plane. Thus the polynomial is hyperbolic and
the source's **closed** cone on each plane is exactly \(L_\theta\succeq0\),
including singular matrices. This is an exact polynomial argument, not an
approximation of the closed boundary.

## 2. An explicit real test basis and its moment domain

Write a test vector in its five ordered four-dimensional blocks as

\[
q=(\Re a,\Im a,\Re b,\Im b,u),\quad
a=A+B w^{-2},\quad b=C+D w^2,
\]

with \(A,B,C,D\in\mathbb C^4\) arbitrary and \(u\in\mathbb R^4\).
Here the coefficient named \(D\) is an arbitrary complex Fourier
coefficient; it is unrelated to the real symmetric matching matrix also
traditionally named \(D\) below.

An explicit real basis consists of the following nine four-element groups.
For each row take \(v=e_k\), \(k=1,2,3,4\); all undisplayed blocks are zero.

| Group | First block | Second block | Third block | Fourth block | Last block |
|---|---|---|---|---|---|
| \(\Re A\) | \(v\) | 0 | 0 | 0 | 0 |
| \(\Im A\) | 0 | \(v\) | 0 | 0 | 0 |
| \(\Re B\) | \(\cos2\theta\,v\) | \(-\sin2\theta\,v\) | 0 | 0 | 0 |
| \(\Im B\) | \(\sin2\theta\,v\) | \(\cos2\theta\,v\) | 0 | 0 | 0 |
| \(\Re C\) | 0 | 0 | \(v\) | 0 | 0 |
| \(\Im C\) | 0 | 0 | 0 | \(v\) | 0 |
| \(\Re D\) | 0 | 0 | \(\cos2\theta\,v\) | \(\sin2\theta\,v\) | 0 |
| \(\Im D\) | 0 | 0 | \(-\sin2\theta\,v\) | \(\cos2\theta\,v\) | 0 |
| \(u\) | 0 | 0 | 0 | 0 | \(v\) |

The complex Fourier coefficients are unique, so these 36 real functions are
independent. Their paired top-block expression is

\[
(\Re a)^T\mathcal X\Re a+(\Im a)^T\mathcal X\Im a
=\Re(a^*\mathcal X a).
\]

Expanding \(a=A+Bw^{-2}\) produces only frequencies \(0,\pm2\): the
\(B^*\mathcal X B\) term is constant, not frequency four. The same
calculation applies to \(b\). The \(\mathcal Z\) term is constant.

Column one of \(U\) has frequencies \(-1,-3\); column two has
frequencies \(+1,+3\). Thus \(a^*U_1\) and \(b^*U_2\), followed by
real parts, use only \(\pm1,\pm3\). The constant columns of \(V\) yield
only \(0,\pm2\). Consequently every quadratic test lies in

\[
\begin{aligned}
\mathcal E_{47}={}&\operatorname{span}\{\mathcal X_{ij},
\mathcal X_{ij}\cos2\theta,\mathcal X_{ij}\sin2\theta\}\;\oplus\;
\operatorname{span}\{\mathcal Z_{ij}\}\\
&\oplus\operatorname{span}\{\eta\cos\theta,\eta\sin\theta,
\eta\cos3\theta,\eta\sin3\theta\}\\
&\oplus\operatorname{span}\{\xi,\xi\cos2\theta,\xi\sin2\theta\},
\qquad 1\leq i\leq j\leq4.
\end{aligned}
\]

Distinct slice-coordinate functions and distinct displayed real harmonics
are independent. Its dimension is \(30+10+4+3=47\).

For the displayed real basis \(q_1,\ldots,q_{36}\), set
\(H_{ij}(\Lambda)=\Lambda(q_i^TL_\theta q_j)\). This is symmetric since
\(L_\theta\) is symmetric. Off-diagonal entries obey the polarization rule

\[
q_i^TLq_j=\tfrac12((q_i+q_j)^TL(q_i+q_j)-q_i^TLq_i-q_j^TLq_j).
\]

In particular, no extra factor of two belongs in \(H_{ij}\); the two
off-diagonal occurrences supply it in \(t^THt\). Polarization also proves
every such entry is in \(\mathcal E_{47}\). The 23 outputs are precisely
the distinct coordinates \(\Lambda(\mathcal X_{ij})\),
\(\Lambda(\mathcal Z_{ij})\), \(\Lambda(\eta\cos\theta)\),
\(\Lambda(\eta\sin\theta)\), \(\Lambda(\xi)\). This leaves 24 free
functional coordinates before any redundancy audit.

## 3. The matching factors always lie in this test space

Fix an arbitrary real symmetric \(S\), with no sign assumption, and
\(z=\tau(\cos\alpha,\sin\alpha,0)+\rho(0,0,1)\), \(\tau>0\).
Let \(w_0=e^{i\alpha}\), \(R_0=\tau U(\alpha)+\rho V\), \(N=SR_0\),
\(v=N_3-iN_2\), \(g=v_++\overline{v_-}\),
\(h=v_+-\overline{v_-}\), and \(\delta=|g|^2+|h|^2\).

For \(\delta>0\), define

\[
a=\frac{|g|^2+|h|^2(w_0/w)^2}{\delta},\qquad
\overline b=\frac{gh(1-(w_0/w)^2)}{\delta},\qquad
M=\begin{pmatrix}a&b\\-\overline b&\overline a\end{pmatrix}.
\]

The checker verifies the linear system for \((a,\overline b)\),
\(MM^*=I_2\), \(M(\alpha)=I_2\), and the two identities

\[
(vM)_++\overline{(vM)_-}=g,\qquad
(vM)_+-\overline{(vM)_-}=h(w_0/w)^2.
\]

These are generic rational identities in unconstrained complex \(g,h\)
and unit-circle Laurent variables. They do not use definiteness of \(S\).

The crucial sharper support statement follows directly from the two matrix
columns:

\[
C_1=(R_0)_1a-(R_0)_2\overline b,
\qquad C_2=(R_0)_1b+(R_0)_2\overline a.
\]

Thus \(C_1\) has only frequencies \(0,-2\), and \(C_2\) only \(0,+2\).
Since \(M\) is unitary, \(\widetilde C\widetilde C^T=Q(z)\).

For completeness, the frequency and derivative identification used by the
source can be checked without an unjustified factor-continuity assertion.
For any complex \(3\times2\) matrix \(T\), put
\(v_T=T_3-iT_2\), \(m_T=T_3+iT_2\), \(n_T=T_1\). Exact expansion gives

\[
\sqrt2\Re\operatorname{tr}(U^*T)=\Re\left[
\frac{w^3}2(v_{T,+}-\overline{v_{T,-}})
+w\left(\frac{m_{T,+}}2-n_{T,+}
-\frac{\overline{m_{T,-}}}2-\overline{n_{T,-}}\right)\right].
\]

For \(T=NM\), the first parenthesis is \(h(w_0/w)^2\). Each term in
the second parenthesis has frequencies \(0,-2\). This leaves only the
first real harmonic. Similarly
\(\sqrt2\langle V,C\rangle_S=\Re g\), a constant.

At \(\alpha\), angular polarization gives the desired values
\(\langle V,R_0\rangle_S=\operatorname{tr}(SQ(k,z))\) and
\(\langle U(\alpha),R_0\rangle_S=\operatorname{tr}(SQ(d(\alpha),z))\).
The checker verifies both for completely generic symmetric \(S\).
Constant \(\widetilde C\widetilde C^T\) implies
\(\langle R_0,C'(\alpha)\rangle_S=0\), and the constant \(V\)-pairing
implies \(\langle V,C'(\alpha)\rangle_S=0\). Since \(\tau>0\),
\(\langle U(\alpha),C'(\alpha)\rangle_S=0\). One can also see the first
stationarity identity from \(M'(\alpha)^*=-M'(\alpha)\) and the
Hermitian matrix \(R_0^*SR_0\); this identity is checked universally for
a generic Hermitian matrix, including indefinite ones.

The remaining derivative is
\(\langle U'(\alpha),R_0\rangle_S
=\operatorname{tr}(SQ(d'(\alpha),z))\), checked symbolically with generic
symmetric \(S\). Value and derivative at one angle determine a first
harmonic. Therefore, for all \(\eta,\xi,\theta\),

\[
\langle F_\theta(\eta,\xi),C(\theta)\rangle_S
=\operatorname{tr}(SQ(Y_\theta(\eta,\xi),z)).
\]

### The separate branch \(\delta=0\)

Here \(g=h=0\), hence \(v=0\). Choose \(M=I_2\) and \(C=R_0\).
This is an everywhere defined constant factor with the required restricted
column supports. The cubic term in the displayed expansion vanishes because
\(v=0\); the rest is a first harmonic. The \(V\)-pairing is constant.
The same value-and-derivative identification applies (and \(C'=0\)).
There is no division by zero and no limiting argument between the two
branches. The checker treats this branch separately and includes the
nonzero indefinite example \(S=\operatorname{diag}(1,-1,0)\),
\(z=(1,0,0)\).

Finally let \(D\in\operatorname{Sym}^4\) be any real symmetric matrix and
\(u\in\mathbb R^4\), and apply matching with \(S=J_uDJ_u^T\). The
matching test's top complex columns are
\(-DJ_u^TC_1\) and \(-DJ_u^TC_2\). Multiplication by this constant real
matrix preserves both one-sided supports and the realification order.
Thus **every** matching test used below belongs to the displayed
36-dimensional space, for arbitrary \(u,D,z\) with nonaxis \(z\).

## 4. Both implications, including every singular point

### A cone point gives a finite functional witness

For \((X,Z,y)\in K\), choose a slice containing \(y\). The determinant
identity gives \(L_{\theta_0}(X,Z,\eta_0,\xi_0)\succeq0\), even when
\(X\) is singular. Evaluation of the 47 functions at this one point is a
finite, real functional with the prescribed outputs. If \(P(\theta_0)\)
has the real test basis as its 36 columns, its matrix is exactly
\(H=P(\theta_0)^TL_{\theta_0}P(\theta_0)\succeq0\). No closure of an
image is used to obtain this witness.

### Every feasible functional gives scalar certificates

Suppose the outputs are prescribed and \(H(\Lambda)\succeq0\).
Constant tests on any one top block belong to the restricted space, so
\(v^TXv\geq0\) for every \(v\), hence \(X\succeq0\).

Fix arbitrary real \(u\), symmetric real \(D\), and nonaxis \(z\).
Use the matching factor just proved and let \(b_j\) be its four ordered
real columns. The test with top blocks \(-DJ_u^Tb_j\) and bottom block
\(u\) belongs to the restricted space. Before applying \(\Lambda\),
its two identities are

\[
\sum_j b_j^TJ_uD\mathcal XDJ_u^Tb_j
=\operatorname{tr}(J_uD\mathcal XDJ_u^TQ(z))
\]

and

\[
-2\sum_j b_j^TJ_uDJ_u^Tr_j
=-2\operatorname{tr}(J_uDJ_u^TQ(Y_\theta(\eta,\xi),z)).
\]

The first uses only constant \(\mathcal X\)-coordinates, and the second
only the three output functions for \(y\). Thus positivity gives

\[
u^TZu+\operatorname{tr}(J_uDXDJ_u^TQ(z))
-2\operatorname{tr}(J_uDJ_u^TQ(y,z))\geq0.\tag{*}
\]

For fixed \(u,D,X,Z,y\), this is a polynomial of degree at most two in
\(z\). It is nonnegative on the dense set \((z_1,z_2)\ne(0,0)\), so
continuity of this **scalar polynomial** gives (*) for every axis point
as well. Nothing asserts that matching matrices \(C\) converge or are
continuous as \(z\) approaches the axis. No feasible functional, no cone
point, and no semidefinite image is being replaced by a limit.

### All four range conditions

Fix a slice representing the actual output \(y\), and let
\(B_j=B(r_j)\) for its four ordered real factor columns. For any
\(v\in\ker X\), put \(D=t vv^T\) and \(z=y\) in (*), for every real
\(t\). Since \(DXD=0\), one obtains

\[
u^TZu-2t(J_uv)^TQ(y)(J_uv)\geq0\qquad(t\in\mathbb R).
\]

Its coefficient of \(t\) must vanish. The exact angular Gram identity
then gives

\[
0=(J_uv)^TQ(y)(J_uv)=\sum_{j=1}^4(v^TB_ju)^2.
\]

Every summand is a square of a real number. Each is zero. Because \(u\)
and \(v\in\ker X\) were arbitrary, this proves separately

\[
\operatorname{range}B_j\subseteq(\ker X)^\perp
=\operatorname{range}X,\qquad j=1,2,3,4.
\]

These conditions cannot be dropped on the singular boundary.

### Pseudoinverse inequality and completed squares

Define the symmetric \(X^\dagger\) to invert positive eigenvalues of
\(X\) and vanish on its kernel. Set \(D=X^\dagger\), \(z=y\) in (*).
The exact relation \(X^\dagger XX^\dagger=X^\dagger\) gives

\[
Z\succeq\Phi_y(X^\dagger)
=\sum_{j=1}^4B_j^TX^\dagger B_j.
\]

The four range conditions give \(XX^\dagger B_j=B_j\). Consequently,
for completely arbitrary \(a_1,\ldots,a_4,u\in\mathbb R^4\),

\[
\begin{aligned}
&\sum_j a_j^TXa_j+2\sum_j a_j^TB_ju+u^TZu\\
&=\sum_j(a_j+X^\dagger B_ju)^TX(a_j+X^\dagger B_ju)
+u^T\left(Z-\sum_jB_j^TX^\dagger B_j\right)u\geq0.
\end{aligned}
\]

This proves the actual slice pencil positive semidefinite and hence the
original degree-20 closed-cone membership. It applies to every singular
\(X\), including \(X=0\), and to every output \(y\), including the axis
and zero. The functional was arbitrary throughout; it was never assumed
to be integration against a measure.

## 5. One exact common-kernel certificate

For the natural 36-element basis above, prescribe constant \(X=I_4\),
constant \(Z=I_4\), and zero values for every other functional coordinate.
The resulting matrix is exactly \(I_{36}\). Equivalently, it is the
average of evaluations with \(\mathcal X=\mathcal Z=I_4\),
\(\eta=\xi=0\), at angles \(0\) and \(\pi/2\). The two paired Fourier
blocks have opposite signs at these angles, canceling their cross terms;
every diagonal basis norm is one.

Thus an explicit sum of the eight constant diagonal output coefficient
matrices is \(I_{36}\). Their full coefficient family has common kernel
\(\{0\}\). This certificate is specific to this pencil. Auxiliary-map
injectivity and its exact rank are covered by the primary redundancy audit,
not asserted here without its coefficient computation.

## 6. Independent exact checker and limitations

`check_matching.py` uses exact SymPy arithmetic with formal complex
conjugation on unit-circle Laurent variables. It does not import the main
coefficient generator. Its generic identities cover arbitrary real symmetric
\(S\) and \(W\); numerical matrices do not stand in for these universal
claims. It additionally exercises a nonzero indefinite \(D\), a nonzero
indefinite \(S\) in the zero-denominator branch, and exact rank-three
boundary witnesses on and off the axis. The example witnesses check all four
range relations and exhibit a rectangular real Gram factor of their
size-20 slice matrices.

Replay with Python and SymPy 1.14.0:

```sh
python audits/check_matching.py
```

The recorded `check_matching.log` contains 42 successful exact assertions.
The universal symbolic calculations and the mathematical argument above
are complementary: the script is not a formal proof-assistant certificate,
and the finite example checks are supplementary regressions. No numerical
tolerance, genericity assumption on \(X\), or measure representation is
part of the verified conclusion.
