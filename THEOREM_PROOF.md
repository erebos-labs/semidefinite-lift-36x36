# An exact 36 × 36 lift with 24 auxiliary variables

Version **1.0.0** · 12 October 2026

This note compresses and independently checks the angular-certificate construction
in OpenAI, *An Exact Semidefinite Lift of a Nonspectrahedral Hyperbolicity Cone*
(5 October 2026). The mathematical source is pinned to
`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb` in
[openai/math](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb).
All seven requested source files are preserved unchanged under `source/`, with
their URLs and SHA-256 hashes in `source/manifest.json`.

**Result.** The proposed 36/24 construction works exactly. Its 47-dimensional
moment coefficient map is injective. The full coefficient family has no common
matrix kernel, and the 24-dimensional auxiliary coefficient map is injective.
Consequently the specified bounded redundancy audit removes no directions.
These are properties of this pencil, not lower bounds on all possible lifts.

The proof below is self-contained in elementary real/complex linear algebra,
polynomial identities, and continuity of scalar polynomials. The accompanying
programs verify exact finite algebraic identities and coefficient/rank data;
they are not a Lean or other proof-assistant formalization of the theorem.

## 1. The exact cone and the pinned baseline

Write `Sym⁴` for real symmetric 4 × 4 matrices. For
\(y=(y_1,y_2,y_3)\), define

\[
Q(y)=\begin{pmatrix}
y_1^2+y_2^2&-y_1y_2&-y_1y_3\\
-y_1y_2&y_2^2+y_3^2&-y_2y_3\\
-y_1y_3&-y_2y_3&y_3^2+y_1^2
\end{pmatrix}.
\]

For \(W=\begin{pmatrix}a&r^T\\r&T\end{pmatrix}\in\operatorname{Sym}^4\), put

\[
\Phi_y(W)=\begin{pmatrix}
\operatorname{tr}(Q(y)T)&-r^TQ(y)\\
-Q(y)r&aQ(y)
\end{pmatrix}.
\]

The polynomial and direction are exactly

\[
p(X,Z,y)=\det\bigl((\det X)Z-\Phi_y(\operatorname{adj}X)\bigr),
\qquad e=(I_4,I_4,0).
\]

Each entry of the matrix inside the determinant is homogeneous of degree five:
\((\det X)Z\) has degree \(4+1\), and the other term has degree \(2+3\).
Because \(p(e)=1\), this is a nonzero homogeneous polynomial of degree **20**.
Throughout this note,

\[
K=\{x=(X,Z,y):\text{every root of }t\mapsto p(te-x)
\text{ is nonnegative}\}.
\tag{1}
\]

Zero roots, their multiplicities, and singular \(X\) are included. No factor of
\(p\) is removed, and no degree-16 polynomial is substituted.

The ambient dimension is \(10+10+3=23\). Output coordinates are the upper
triangular entries of \(X\), then those of \(Z\), each in lexicographic order,
then \(y_1,y_2,y_3\).

The pinned `introduction.tex` states a homogeneous 100 × 100 real symmetric
pencil with 307 auxiliary variables for this entire closed cone.
The count in `lift.tex` is \(20\dim\mathcal T_2=20\cdot5=100\) test
coordinates and \(22\dim\mathcal T_7=22\cdot15=330\) functional coordinates,
of which 23 are prescribed outputs. Thus \(330-23=307\). The proof in the
source supplies the singular range conditions and the pseudoinverse completion.
Sections 5–8 below reprove the ingredients needed for the smaller restriction.

The separate family-095 manuscript *Hyperbolicity Cones Without Semidefinite
Lifts* concerns different examples. Its nonliftability assertion is not an
assertion about the fixed cone (1), and is not an input to this proof.

## 2. The explicit pencil and theorem

Use auxiliary matrices \(X_c,X_s\in\operatorname{Sym}^4\) and four real
auxiliary scalars \(h_c,h_s,k_c,k_s\). These comprise \(10+10+4=24\)
independent real coordinates.

For \(\varepsilon\in\{1,-1\}\), define the 16 × 16 matrix in 4 × 4 blocks

\[
T_\varepsilon(X,X_c,X_s)=
\begin{pmatrix}
X&0&X_c&\varepsilon X_s\\
0&X&-\varepsilon X_s&X_c\\
X_c&-\varepsilon X_s&X&0\\
\varepsilon X_s&X_c&0&X
\end{pmatrix}.
\tag{2}
\]

For a column \(r\in\mathbb R^3\), retain the source's sign convention

\[
B(r)=\begin{pmatrix}0&r^T\\-r&0_{3\times3}\end{pmatrix}.
\tag{3}
\]

Define eight columns \(\ell_j\in\mathbb R^3\) by the following table. The
listed triple is **\(\sqrt2\ell_j\)**, so every row has the same factor
\(1/\sqrt2\).

| \(j\) | \(\sqrt2\ell_j\) |
|---|---|
| 1 | \(\bigl(-y_1,(h_s-y_2)/2,(h_c+y_1)/2+y_3\bigr)^T\) |
| 2 | \(\bigl(y_2,(h_c-y_1)/2+y_3,-(h_s+y_2)/2\bigr)^T\) |
| 3 | \(\bigl(-y_1,y_2-k_s,y_1+k_c\bigr)^T\) |
| 4 | \(\bigl(-y_2,k_c,k_s\bigr)^T\) |
| 5 | \(\bigl(-y_1,(h_s-y_2)/2,-(h_c+y_1)/2+y_3\bigr)^T\) |
| 6 | \(\bigl(-y_2,(y_1-h_c)/2+y_3,-(h_s+y_2)/2\bigr)^T\) |
| 7 | \(\bigl(-y_1,y_2+k_s,-y_1+k_c\bigr)^T\) |
| 8 | \(\bigl(y_2,k_c,-k_s\bigr)^T\) |

Stack these skew matrices as

\[
K_a=\begin{pmatrix}B(\ell_1)\\B(\ell_2)\\B(\ell_3)\\B(\ell_4)\end{pmatrix},
\qquad
K_b=\begin{pmatrix}B(\ell_5)\\B(\ell_6)\\B(\ell_7)\\B(\ell_8)\end{pmatrix}.
\]

The complete homogeneous pencil is

\[
\mathcal H(X,Z,y;X_c,X_s,h_c,h_s,k_c,k_s)=
\begin{pmatrix}
T_+(X,X_c,X_s)&0&K_a\\
0&T_-(X,X_c,X_s)&K_b\\
K_a^T&K_b^T&Z
\end{pmatrix}\in\operatorname{Sym}^{36}.
\tag{4}
\]

Despite the skewness of each \(B(\ell_j)\), the whole pencil is symmetric:
the lower off-diagonal blocks in (4) are the **transposes** of the upper ones.

### Theorem

For every \((X,Z,y)\in\operatorname{Sym}^4\times\operatorname{Sym}^4
\times\mathbb R^3\),

\[
(X,Z,y)\in K
\quad\Longleftrightarrow\quad
\exists\,X_c,X_s,h_c,h_s,k_c,k_s\quad
\mathcal H(X,Z,y;X_c,X_s,h_c,h_s,k_c,k_s)\succeq0.
\tag{5}
\]

Equivalently there are explicitly specified
\(A_1,\ldots,A_{23},B_1,\ldots,B_{24}\in
\operatorname{Sym}^{36}(\mathbb Q(\sqrt2))\) with

\[
K=\left\{x:\exists v\in\mathbb R^{24},\quad
\sum_{\nu=1}^{23}x_\nu A_\nu+\sum_{j=1}^{24}v_jB_j\succeq0\right\}.
\]

There is no constant matrix term. Equality holds as written, with no closure
of the projection. All coefficients of (4) belong to
\(\{0,\pm1,\pm\sqrt2/2,\pm\sqrt2/4\}\).

## 3. Ordered real test basis and moment coordinates

Let \(w=e^{i\theta}\). Use the source realification order

\[
q(\theta)=(\operatorname{Re}a,\operatorname{Im}a,
             \operatorname{Re}b,\operatorname{Im}b,u),
\quad a=A+B w^{-2},\quad b=C+D w^2,
\tag{6}
\]

where \(A,B,C,D\in\mathbb C^4\) and \(u\in\mathbb R^4\).
Here \(D\) is a complex coefficient vector; the symmetric matrix named \(D\)
in Sections 7–8 is a different, explicitly indicated object.

For \(i=0,1,2,3\), let \(e_i\) be the corresponding standard basis vector of
\(\mathbb R^4\), and put \(c_2=\cos2\theta,s_2=\sin2\theta\).
The following table is the complete real basis \(q_1,\ldots,q_{36}\): within
each row, take \(i=0,1,2,3\) in order.

| Basis numbers (one based) | Coefficient group | Five 4-vector blocks |
|---|---|---|
| 1–4 | \(\operatorname{Re}A\) | \((e_i,0,0,0,0)\) |
| 5–8 | \(\operatorname{Im}A\) | \((0,e_i,0,0,0)\) |
| 9–12 | \(\operatorname{Re}B\) | \((c_2e_i,-s_2e_i,0,0,0)\) |
| 13–16 | \(\operatorname{Im}B\) | \((s_2e_i,c_2e_i,0,0,0)\) |
| 17–20 | \(\operatorname{Re}C\) | \((0,0,e_i,0,0)\) |
| 21–24 | \(\operatorname{Im}C\) | \((0,0,0,e_i,0)\) |
| 25–28 | \(\operatorname{Re}D\) | \((0,0,c_2e_i,s_2e_i,0)\) |
| 29–32 | \(\operatorname{Im}D\) | \((0,0,-s_2e_i,c_2e_i,0)\) |
| 33–36 | \(u\) | \((0,0,0,0,e_i)\) |

The coefficients of \(1,w^{-2}\) in \(a\), those of \(1,w^2\) in \(b\),
and the constant \(u\) are uniquely determined. Therefore these 36 vectors
are independent as vector-valued functions. Denote their real span by
\(\mathcal W\), and let \(P(\theta)\) be the 20 × 36 matrix with these
columns.

Use formal slice variables \(\mathcal X,\mathcal Z,\eta,\xi\). Define
the real vector space of scalar functions

\[
\begin{aligned}
\mathcal E_{\rm small}={}&
\operatorname{span}\{\mathcal X_{ij},\mathcal X_{ij}\cos2\theta,
\mathcal X_{ij}\sin2\theta:i\le j\}\\
&\oplus\operatorname{span}\{\mathcal Z_{ij}:i\le j\}\\
&\oplus\operatorname{span}\{\eta\cos\theta,\eta\sin\theta,
\eta\cos3\theta,\eta\sin3\theta\}\\
&\oplus\operatorname{span}\{\xi,\xi\cos2\theta,\xi\sin2\theta\}.
\end{aligned}
\tag{7}
\]

Independent slice coordinate functions and independent real trigonometric
functions make this a direct sum with dimension
\(10\cdot3+10+4+3=47\).

A real linear functional \(\Lambda\) on (7) is ordered as follows:

1. For each upper pair \((i,j)\) in lexicographic order, list the moments
   \(\Lambda(\mathcal X_{ij}),\Lambda(\mathcal X_{ij}\cos2\theta),
   \Lambda(\mathcal X_{ij}\sin2\theta)\).
2. List the ten \(\Lambda(\mathcal Z_{ij})\) in the same pair order.
3. List \(\Lambda(\eta\cos\theta),\Lambda(\eta\sin\theta),
   \Lambda(\eta\cos3\theta),\Lambda(\eta\sin3\theta)\).
4. List \(\Lambda(\xi),\Lambda(\xi\cos2\theta),\Lambda(\xi\sin2\theta)\).

Prescribe

\[
\Lambda(\mathcal X)=X,\quad \Lambda(\mathcal Z)=Z,\quad
(\Lambda(\eta\cos\theta),\Lambda(\eta\sin\theta),\Lambda(\xi))=y.
\tag{8}
\]

These are 23 distinct basis coordinates. Their selector matrix has an identity
23 × 23 minor, so independence does not rely on an informal dimension count.
The remaining coordinates are exactly

\[
X_c=\Lambda(\mathcal X\cos2\theta),\quad
X_s=\Lambda(\mathcal X\sin2\theta),\quad
h_c=\Lambda(\eta\cos3\theta),\ h_s=\Lambda(\eta\sin3\theta),\quad
k_c=\Lambda(\xi\cos2\theta),\ k_s=\Lambda(\xi\sin2\theta).
\]

In the executable order, auxiliary \(X_c,X_s\) entries are interleaved per
upper triangular pair, followed by \(h_c,h_s,k_c,k_s\). There are 24.

## 4. Coefficient generation and the frequency audit

The source slice pencil \(L_\theta\) is defined in Section 5 below. Set

\[
H(\Lambda)_{ij}=\Lambda\bigl(q_i^T L_\theta q_j\bigr),
\qquad 1\le i,j\le36.
\tag{9}
\]

For clarity, the off-diagonal entry is the bilinear polarization

\[
q_i^TL_\theta q_j=
\tfrac12\bigl((q_i+q_j)^TL_\theta(q_i+q_j)
-q_i^TL_\theta q_i-q_j^TL_\theta q_j\bigr).
\]

Thus there is no extra factor of two in a stored entry; the factor two for
cross terms appears when evaluating \(t^THt\).

The paired top-block energy is

\[
(\operatorname{Re}a)^T\mathcal X\operatorname{Re}a+
(\operatorname{Im}a)^T\mathcal X\operatorname{Im}a
=\operatorname{Re}(a^*\mathcal X a),
\]

and similarly for \(b\). Inserting (6) leaves only frequencies \(0,\pm2\):
the pure \(B\) contribution has \(w^2w^{-2}=1\), and the apparent fourth
harmonics cancel between the real and imaginary blocks. Polarizing this
identity gives all top entries, not only the diagonal tests. It produces the
two blocks (2). The bottom vector is constant, so \(\mathcal Z\) appears
only at frequency zero.

Write \(F=(f_1,f_2)=\eta U+\xi V_0\), as in Section 5. The eight top/bottom
cross blocks, before applying \(\Lambda\), are \(B\) applied to

\[
\operatorname{Re}f_1,\ \operatorname{Im}f_1,\
\operatorname{Re}(w^2f_1),\ \operatorname{Im}(w^2f_1),\
\operatorname{Re}f_2,\ \operatorname{Im}f_2,\
\operatorname{Re}(w^{-2}f_2),\ \operatorname{Im}(w^{-2}f_2).
\tag{10}
\]

In particular,

\[
\sqrt2\,w^2f_1=\eta(-w,s,c)^T+\xi(0,iw^2,w^2)^T,
\]

\[
\sqrt2\,w^{-2}f_2=\eta(-w^{-1},s,-c)^T+\xi(0,iw^{-2},w^{-2})^T.
\]

The unshifted factors have \(\eta\)-frequencies \(\pm1,\pm3\) and
constant \(\xi\)-terms. The shifted factors have first harmonics for
\(\eta\) and second harmonics for \(\xi\). There are no fifth harmonics.
Expanding (10), using
\(s\cos2\theta=(\sin3\theta-\sin\theta)/2\) and the other elementary
product identities, gives exactly the eight columns in the table in Section 2.

This proves that every product in (9) belongs to (7) and that (9) equals (4).
It is stronger than frequency counting: the complete explicit coefficient
rule is (2)–(4), and `audits/independent_coefficients.py` separately multiplies
all source Laurent expressions in (9), including polarization and the stated
realification order. The resulting 47 matrices agree entry by entry with
`generator.py` and `coefficients.json`.

Finally, because \(q=\sum_i t_iq_i\),

\[
H(\Lambda)\succeq0\quad\Longleftrightarrow\quad
\Lambda(q^TL_\theta q)\ge0\ \text{for all }q\in\mathcal W.
\tag{11}
\]

No representing measure is part of the definition of \(\Lambda\).

## 5. Exact angular slices, including their singular points

Put \(c=\cos\theta,s=\sin\theta,w=e^{i\theta}\), and define

\[
d(\theta)=(c,s,0),\quad k=(0,0,1),\quad
Y_\theta(\eta,\xi)=\eta d(\theta)+\xi k.
\]

These planes cover \(\mathbb R^3\). Retain the source factors

\[
U(\theta)=\frac1{\sqrt2}
\begin{pmatrix}-w^{-1}&-w\\sw^{-2}&sw^2\\cw^{-2}&-cw^2\end{pmatrix},
\qquad
V_0=\frac1{\sqrt2}\begin{pmatrix}0&0\\i&i\\1&1\end{pmatrix},
\qquad F_\theta=\eta U(\theta)+\xi V_0.
\tag{12}
\]

For \(F=(f_1,f_2)\in\mathbb C^{3\times2}\), realification means

\[
\widetilde F=(\operatorname{Re}f_1,\operatorname{Im}f_1,
               \operatorname{Re}f_2,\operatorname{Im}f_2).
\]

It satisfies \(\widetilde F\widetilde F^T=\operatorname{Re}(FF^*)\).
Multiplication of (12) gives

\[
\operatorname{Re}(UU^*)=
\begin{pmatrix}1&-sc&0\\-sc&s^2&0\\0&0&c^2\end{pmatrix},\quad
\operatorname{Re}(V_0V_0^*)=\operatorname{diag}(0,1,1),
\]

\[
\operatorname{Re}(UV_0^*+V_0U^*)=
\begin{pmatrix}0&0&-c\\0&0&-s\\-c&-s&0\end{pmatrix}.
\]

Consequently

\[
\widetilde F_\theta\widetilde F_\theta^T
=Q(Y_\theta(\eta,\xi)).
\tag{13}
\]

In particular \(Q(y)\succeq0\) for every real \(y\).
For \(u=(u_0,u')\in\mathbb R\times\mathbb R^3\), set

\[
J_u=\begin{pmatrix}u'&-u_0I_3\end{pmatrix}.
\]

Direct multiplication proves \(B(r)u=J_u^Tr\). If \(r_1,\ldots,r_4\)
are the ordered columns of \(\widetilde F_\theta\), then for every symmetric
\(W\), including indefinite \(W\),

\[
\sum_{j=1}^4 B(r_j)^TWB(r_j)=\Phi_{Y_\theta}(W),\qquad
u^T\Phi_y(W)u=\operatorname{tr}(J_uWJ_u^TQ(y)).
\tag{14}
\]

For example, the single-column summand for
\(W=\begin{pmatrix}a&v^T\\v&T\end{pmatrix}\) has blocks
\(r^TTr,-r^Tvr^T,-rr^Tv,arr^T\); summation uses (13).

Let \(\mathcal B\) be the 16 × 4 vertical stack of these four \(B(r_j)\).
The original 20 × 20 slice pencil is

\[
L_\theta(\mathcal X,\mathcal Z,\eta,\xi)=
\begin{pmatrix}I_4\otimes\mathcal X&\mathcal B\\
\mathcal B^T&\mathcal Z\end{pmatrix}.
\tag{15}
\]

Here is an exact polynomial proof of its determinant identity which does not
assume that an evaluated \(\mathcal X\) is invertible. In the polynomial ring
of slice variables, write \(d=\det\mathcal X\),
\(A=I_4\otimes\mathcal X\), and
\(\widehat J=I_4\otimes\operatorname{adj}\mathcal X\). Then
\(A\widehat J=dI_{16}\). Multiply (15) on the right by

\[
R=\begin{pmatrix}dI_{16}&-\widehat J\mathcal B\\0&dI_4\end{pmatrix}.
\]

Using (14) gives the polynomial matrix identity

\[
L_\theta R=
\begin{pmatrix}
dA&0\\
d\mathcal B^T&d\mathcal Z-\Phi_{Y_\theta}(\operatorname{adj}\mathcal X)
\end{pmatrix}.
\tag{16}
\]

Now \(\det R=d^{20}\) and
\(\det(dA)=d^{16}\det A=d^{20}\). Taking determinants in (16) yields
\(d^{20}\det L_\theta=d^{20}p(\mathcal X,\mathcal Z,Y_\theta)\).
Since \(d\) is a nonzero polynomial, cancellation in the polynomial ring gives

\[
\det L_\theta=p(\mathcal X,\mathcal Z,Y_\theta)
\tag{17}
\]

as a polynomial identity. Evaluating this identity is valid even when
\(\det\mathcal X=0\). This cancellation proves (17); it does not remove a
factor from the defining polynomial \(p\).

Moreover \(L_\theta(I_4,I_4,0,0)=I_{20}\). Linearity of the pencil and of
\(Y_\theta\) therefore gives, at every point of the slice,

\[
p\bigl(te-(\mathcal X,\mathcal Z,Y_\theta(\eta,\xi))\bigr)
=\det\bigl(tI_{20}-L_\theta(\mathcal X,\mathcal Z,\eta,\xi)\bigr).
\tag{18}
\]

The right side is the characteristic polynomial of a real symmetric matrix.
Its roots are its twenty real eigenvalues with multiplicities. Thus \(p\) is
hyperbolic at \(e\), and definition (1) is equivalent to

\[
(\mathcal X,\mathcal Z,Y_\theta(\eta,\xi))\in K
\quad\Longleftrightarrow\quad L_\theta\succeq0.
\tag{19}
\]

This proves slice equality directly, including all zero roots and all singular
\(\mathcal X\). There is no limiting positive-definite approximation in this
argument.

## 6. Matching certificates and their restricted support

Define symmetric bilinear polarization

\[
Q(y,z)=\tfrac12\bigl(Q(y+z)-Q(y)-Q(z)\bigr).
\]

For real symmetric \(S\in\operatorname{Sym}^3\), use the real bilinear
pairing \(\langle F,G\rangle_S=\operatorname{Re}\operatorname{tr}(F^*SG)\).
The argument below imposes no sign restriction on \(S\).

**Matching claim.** If \((z_1,z_2)\ne(0,0)\), there is
\(\mathscr C(\theta)=(\mathscr C_1,\mathscr C_2)\) with

\[
\mathscr C_1\in\operatorname{span}_{\mathbb C}\{1,w^{-2}\}^3,
\qquad
\mathscr C_2\in\operatorname{span}_{\mathbb C}\{1,w^2\}^3,
\tag{20}
\]

such that, for every \(\theta,\eta,\xi\),

\[
\widetilde{\mathscr C}\widetilde{\mathscr C}^{T}=Q(z),\qquad
\langle F_\theta,\mathscr C\rangle_S
=\operatorname{tr}\bigl(SQ(Y_\theta,z)\bigr).
\tag{21}
\]

Here is the full construction. Write
\(z=\tau d(\alpha)+\rho k\), with \(\tau>0\), and set

\[
w_0=e^{i\alpha},\quad R_0=\tau U(\alpha)+\rho V_0,\quad N=SR_0,
\quad v=N_3-iN_2.
\]

Rows of \(N\) are numbered 1,2,3; subscripts \(+,-\) below mean the two
components of a row. Put

\[
g=v_++\overline{v_-},\quad h=v_+-\overline{v_-},\quad
\kappa=(w_0/w)^2,\quad\delta=|g|^2+|h|^2.
\]

If \(\delta>0\), define

\[
a_0(\theta)=\frac{|g|^2+|h|^2\kappa}{\delta},\qquad
\overline{b_0(\theta)}=\frac{gh(1-\kappa)}{\delta},\qquad
M=\begin{pmatrix}a_0&b_0\\-\overline b_0&\overline a_0\end{pmatrix}.
\tag{22}
\]

These expressions solve

\[
\begin{pmatrix}g&\overline h\\h&-\overline g\end{pmatrix}
\begin{pmatrix}a_0\\\overline b_0\end{pmatrix}
=\begin{pmatrix}g\\h\kappa\end{pmatrix}.
\]

The coefficient matrix has adjoint-times-itself equal to \(\delta I_2\).
The right side has squared norm \(\delta\), so
\(|a_0|^2+|b_0|^2=1\); hence \(M\) is unitary. At \(\theta=\alpha\),
\(M=I_2\).

If \(\delta=0\), then \(g=h=0\), equivalently \(v=0\). Define \(M=I_2\)
directly. This is a separate valid branch, with no division and no limiting
choice of factors.

In both branches put \(\mathscr C=R_0M\). Formula (22) shows that
\(a_0,\overline b_0\) use only \(0,-2\), whereas their conjugates use
only \(0,+2\). Thus

\[
\mathscr C_1=(R_0)_1a_0-(R_0)_2\overline b_0,\qquad
\mathscr C_2=(R_0)_1b_0+(R_0)_2\overline a_0
\]

proves (20), including the constant \(M=I\) branch. Unitarity and (13) give
the first identity in (21).

To prove the second, write \(r^M=rM\) for any row. The linear system gives

\[
v^M_++\overline{v^M_-}=g,\qquad
v^M_+-\overline{v^M_-}=h\kappa.
\tag{23}
\]

These also hold when \(v=0,M=I\). Directly from \(V_0\),
\(\sqrt2\langle V_0,\mathscr C\rangle_S=\operatorname{Re}g\), a constant.
Writing \(m=N_3+iN_2,n=N_1\), expansion of \(U\) gives

\[
\sqrt2\langle U,\mathscr C\rangle_S=
\operatorname{Re}\left[
\frac{w^3}{2}(v^M_+-\overline{v^M_-})+
w\left(\frac{m^M_+}{2}-n^M_+-
\frac{\overline{m^M_-}}{2}-\overline{n^M_-}\right)\right].
\tag{24}
\]

The first term reduces to frequency one by \(w^3\kappa=w_0^2w\).
In the second, each \((rM)_+\) and \(\overline{(rM)_-}\) has only
frequencies \(0,-2\). Thus the entire real pairing is a real first harmonic.
For \(\delta=0\), the cubic term is zero, and the same conclusion follows.

At \(\alpha\), polarization of (13) yields

\[
\langle V_0,R_0\rangle_S=\operatorname{tr}(SQ(k,z)),\qquad
\langle U(\alpha),R_0\rangle_S=\operatorname{tr}(SQ(d(\alpha),z)).
\]

The constant \(V_0\) pairing is therefore identified. To identify the first
harmonic, differentiate the constant Gram identity for \(\mathscr C\):
\(\langle R_0,\mathscr C'(\alpha)\rangle_S=0\).
Constancy of its \(V_0\) pairing gives
\(\langle V_0,\mathscr C'(\alpha)\rangle_S=0\).
Using \(R_0=\tau U(\alpha)+\rho V_0\) and \(\tau>0\) gives
\(\langle U(\alpha),\mathscr C'(\alpha)\rangle_S=0\).
Differentiating (13) for \(\tau U(\theta)+\rho V_0\) at \(\alpha\) then
gives

\[
\langle U'(\alpha),R_0\rangle_S
=\operatorname{tr}(SQ(d'(\alpha),z)).
\]

Hence the first harmonic in (24) and
\(\operatorname{tr}(SQ(d(\theta),z))\) have the same value and derivative
at \(\alpha\). A real first harmonic is determined by these two numbers,
so they are identical. Linearity in \(\eta,\xi\) proves (21).
All differentiation is of explicitly defined trigonometric functions for
fixed \(S,z\). No positive definiteness of \(S\) was used.

## 7. Cone membership gives a finite feasible witness

Let \((X,Z,y)\in K\). Choose real \(\theta_0,\eta_0,\xi_0\) with
\(y=Y_{\theta_0}(\eta_0,\xi_0)\). For example, if
\(r=\sqrt{y_1^2+y_2^2}>0\), choose
\(\eta_0=r\), \(\cos\theta_0=y_1/r\), \(\sin\theta_0=y_2/r\),
and \(\xi_0=y_3\). If \(r=0\), choose \(\theta_0=0,\eta_0=0,\xi_0=y_3\).

Define \(\Lambda\) by evaluation of every basis function (7) at
\((\theta_0,X,Z,\eta_0,\xi_0)\). These are 47 finite real numbers, with
the prescribed outputs (8). Specifically,

\[
X_c=X\cos2\theta_0,\quad X_s=X\sin2\theta_0,\quad
h_c=\eta_0\cos3\theta_0,\quad h_s=\eta_0\sin3\theta_0,\quad
k_c=\xi_0\cos2\theta_0,\quad k_s=\xi_0\sin2\theta_0.
\tag{25}
\]

By (19), \(L_{\theta_0}(X,Z,\eta_0,\xi_0)\succeq0\). Therefore

\[
H(\Lambda)=P(\theta_0)^T L_{\theta_0}(X,Z,\eta_0,\xi_0)P(\theta_0)
\succeq0.
\]

This proves one direction of (5), with an explicit finite witness at every
point, including singular \(X\). No limiting auxiliary values are needed.

## 8. An arbitrary feasible functional gives cone membership

Suppose a real \(\Lambda\) on (7) has outputs (8) and \(H(\Lambda)\succeq0\).
It need not be an evaluation or a moment functional of a measure.
Constant tests supported in the first 4-vector block are in \(\mathcal W\);
by (11) they imply \(v^TXv\ge0\) for all \(v\), hence \(X\succeq0\).

Fix arbitrary \(u\in\mathbb R^4\), an arbitrary **symmetric matrix**
\(D\in\operatorname{Sym}^4\), and initially
\(z\) with \((z_1,z_2)\ne(0,0)\). Apply Section 6 with
\(S=J_uDJ_u^T\). In particular, \(D\) and \(S\) may be indefinite.
Let \(b_1,\ldots,b_4\) be the source-ordered real columns of
\(\widetilde{\mathscr C}\). Use the test

\[
q(\theta)=(-DJ_u^Tb_1,-DJ_u^Tb_2,-DJ_u^Tb_3,-DJ_u^Tb_4,u).
\tag{26}
\]

The complex first pair equals \(-DJ_u^T\mathscr C_1\), and the second equals
\(-DJ_u^T\mathscr C_2\). The map \(-DJ_u^T\) is constant and real, so (20)
proves that (26) belongs to **the restricted space** \(\mathcal W\).
This containment, not the mere degree bound, is what makes the compression
exact.

Before applying \(\Lambda\), the top contribution is

\[
\sum_j b_j^TJ_uD\mathcal XDJ_u^Tb_j
=\operatorname{tr}(J_uD\mathcal XDJ_u^TQ(z)),
\]

and the cross contribution is

\[
-2\sum_j b_j^TJ_uDJ_u^Tr_j
=-2\langle F_\theta,\mathscr C\rangle_{J_uDJ_u^T}
=-2\operatorname{tr}(J_uDJ_u^TQ(Y_\theta,z)).
\]

The first is angle independent in its \(\mathcal X\) coefficient; the
second is linear in \(Y_\theta\). Both expressions belong to (7).
Applying (8) and (11) yields the scalar inequality

\[
u^TZu+\operatorname{tr}(J_uDXDJ_u^TQ(z))
-2\operatorname{tr}(J_uDJ_u^TQ(y,z))\ge0.
\tag{27}
\]

For fixed \(u,D,X,Z,y\), its left side is a polynomial in \(z\).
If \(z_1=z_2=0\), apply (27) to \(z^{(n)}=(1/n,0,z_3)\) and use
continuity of this scalar polynomial. It follows that (27) holds for **all**
\(z\in\mathbb R^3\). This step makes no assertion about continuity or
limits of the matching factors \(\mathscr C\), and takes no closure of a
projected cone.

Now fix any slice representation \(y=Y_{\theta_0}(\eta_0,\xi_0)\), and
let \(r_j\) be its four real Gram columns, with \(B_j=B(r_j)\).

### All four range conditions

For \(v\in\ker X\), take \(D=t vv^T\), where \(t\in\mathbb R\), and
take \(z=y\) in (27). Then \(DXD=0\), so

\[
u^TZu-2t(J_uv)^TQ(y)(J_uv)\ge0\qquad\text{for every real }t.
\]

The coefficient of \(t\) must vanish. By (13) and \(B_ju=J_u^Tr_j\),

\[
0=(J_uv)^TQ(y)(J_uv)=\sum_{j=1}^4(v^TB_ju)^2.
\]

Each summand is a real square, so every \(v^TB_ju\) is zero, separately
for all four \(j\). As \(u\) and \(v\in\ker X\) are arbitrary,

\[
\operatorname{range}B_j\subseteq(\ker X)^\perp
=\operatorname{range}X\qquad(j=1,2,3,4).
\tag{28}
\]

### The pseudoinverse inequality

Let \(X^\dagger\) be the real symmetric matrix that inverts the positive
eigenvalues of \(X\) and is zero on its kernel. It exists also for \(X=0\)
and satisfies \(X^\dagger XX^\dagger=X^\dagger\).
Put \(D=X^\dagger,z=y\) in (27). Using (14) gives, for every \(u\),

\[
u^T\bigl(Z-\Phi_y(X^\dagger)\bigr)u\ge0,
\quad\text{hence}\quad
Z\succeq\Phi_y(X^\dagger)=\sum_{j=1}^4 B_j^TX^\dagger B_j.
\tag{29}
\]

### Completed squares prove actual slice positivity

For arbitrary \(a_1,\ldots,a_4,u\in\mathbb R^4\), (28) gives
\(XX^\dagger B_j=B_j\). Therefore the exact identity

\[
\begin{aligned}
&\sum_j a_j^TXa_j+2\sum_j a_j^TB_ju+u^TZu\\
&=\sum_j(a_j+X^\dagger B_ju)^TX(a_j+X^\dagger B_ju)
+u^T\left(Z-\sum_jB_j^TX^\dagger B_j\right)u
\end{aligned}
\tag{30}
\]

has a nonnegative right side by \(X\succeq0\) and (29).
Thus \(L_{\theta_0}(X,Z,\eta_0,\xi_0)\succeq0\). Equation (19) proves
\((X,Z,y)\in K\), finishing the converse and the theorem.

Every point allowed by (4) has been shown to satisfy the original root
condition (1). The proof used neither a representing measure nor an assumed
closedness of a spectrahedral projection.

## 9. One bounded redundancy audit

Let \(H_k\), \(0\le k<47\), denote the coefficient matrices in the ordered
moment system of Section 3. These indices are **zero based** in all JSON files.

### Common matrix kernel of the full pencil

Set constant \(X=Z=I_4\) and every other moment to zero in (4). This yields
\(I_{36}\) exactly. More explicitly,

\[
H_0+H_{12}+H_{21}+H_{27}+H_{30}+H_{34}+H_{37}+H_{39}=I_{36}.
\tag{31}
\]

All eight matrices are output coefficient matrices. Thus the common kernel
of the output family, and also of the whole coefficient family, is \(\{0\}\).
There is no common zero subspace available for a congruence compression.

### Auxiliary parameter map

With outputs fixed, the linear map from the 24 auxiliary coordinates to the
upper triangular entries of \(\mathcal H\) has a diagonal 24 × 24 minor.
For each upper pair \((i,j)\), select entries
\((i,8+j)\) for \((X_c)_{ij}\), and \((i,12+j)\) for \((X_s)_{ij}\).
For \(h_c,h_s,k_c,k_s\), respectively select
\((4,34),(0,34),(12,34),(12,35)\).
All these matrix indices are zero based.

In auxiliary order the selected minor is

\[
\operatorname{diag}\left(
1^{\times20},\frac{\sqrt2}{4},\frac{\sqrt2}{4},
\frac{\sqrt2}{2},\frac{\sqrt2}{2}\right),
\quad\det=\frac1{16}\ne0.
\tag{32}
\]

Its rank is therefore 24 and its kernel is zero. No auxiliary direction
can be deleted as a linear redundancy with outputs fixed.

The full moment-to-matrix map likewise has rank 47. Both independent
implementations supply an explicit 47 × 47 minor with determinant
\(\sqrt2/32\); row orderings are recorded with the respective certificates.
Consequently the actual span of all \(q_i^TL_\theta q_j\) is precisely (7),
not a smaller unverified subspace.

The common vector kernel of the **auxiliary coefficient matrices alone** is a
different object from the kernel of the auxiliary parameter map. In this basis
it is \(\operatorname{span}\{e_{33}\}\), the bottom \(u_1\) direction
(zero-based index 33). Indeed the \(X_c=I_4\) auxiliary combination is
invertible on the first 32 coordinates and zero on the last four. On the
remaining four coordinates, the \(k_c\) coefficient forces
\(u_0=u_2=u_3=0\); every auxiliary matrix annihilates \(u_1\).
This direction is not a common kernel of the complete pencil, since output
\(Z\) acts on it. It cannot be removed by the common-kernel criterion.

The audit is confined to these exact kernels and ranks of this fixed
construction. It performs no search over alternative lifts, block
reformulations, or globally minimal matrix sizes.

## 10. Verification record and limits

The replay consists of complementary checks:

| Component | Exact claim checked |
|---|---|
| `generator.py` | Compact formulas (2)–(4) generate 47 homogeneous symmetric coefficient matrices over \(\mathbb Q(\sqrt2)\). |
| `audits/independent_coefficients.py` | Direct source Laurent multiplication in \(\mathbb Q(\sqrt2,i)[w,w^{-1}]\), all 28,512 directed contractions, supports, realification order, and entrywise agreement of all 47 matrices. |
| `verify.py` | Pinned source hashes, saved/generated coefficient equality, the exact kernel/minor certificates, and singular-boundary regression examples. |
| `audits/check_matching.py` | Universal angular, rotation, pairing, symmetric-\(D\) certificate identities; both rotation branches; exact indefinite and singular examples. |
| `audits/check_determinant_identity.py` | A denominator-free polynomial certificate for all 400 entries of (16), adjugate identities, and the four-block \(\Phi\) factorization. |

There are 60,912 entries in the 47 full 36 × 36 matrices; exact symmetry permits
31,302 upper-triangle coefficients, of which 298 are nonzero. The independent
program reconstructs directed lower contractions as well as upper ones. No
numerical tolerance, floating-point rank decision, numerical SDP feasibility
test, generic point, or Monte Carlo sample certifies the theorem.

An exact nonzero-\(y\) singular regression point is

\[
X=\operatorname{diag}(1,0,1,1),\quad
y=(0,0,1),\quad Z=\operatorname{diag}(2,0,1,1).
\]

With \(\theta=0,\eta=0,\xi=1\), the source slice and its lifted witness both
have rank 12 and explicit exact Gram factorizations. This tests boundary
implementation, while the universal boundary proof is Section 8.

To expose the necessity of (28), use instead
\(X=\operatorname{diag}(1,0,0,0)\), the same \(y\), and
\(Z=\Phi_y(X^\dagger)=\operatorname{diag}(0,0,1,1)\).
The pseudoinverse inequality holds with equality, but a source slice principal
minor has determinant \(-1/2\). This point is outside the cone because the
range condition fails. The checker verifies this exactly.

The universal theorem still uses the written linear-algebra argument: the
spectral theorem for real symmetric matrices, polynomial-ring cancellation,
determination of a first harmonic by its value and derivative, scalar
continuity, and the completed-square argument. The scripts are transparent
symbolic verifiers of its finite identities, not a proof-assistant kernel
checking every quantified logical step. No claim is made about numerical
conditioning or practical SDP solver performance.

## 11. Attribution, later changes, and originality

The pinned source already contains the 20 × 20 angular pencil, the unitary
matching construction, and the complete singular-boundary converse. Its
introduction attributes the finite-functional framework to the moment/SOS
literature. This result preserves and reuses those ingredients.

The substantive observation in this compression is that the matching factor
has separate first-column support \(\{0,-2\}\) and second-column support
\(\{0,+2\}\). Keeping the paired real and imaginary coordinates preserves
every matching certificate while cancelling the spurious fourth harmonics
in the top block. Once this restriction is justified, expanding the blocks,
deleting unused moment coordinates, counting 36/47/24, and checking linear
redundancy are finite algebraic work. This is a smaller explicit representation
of this fixed cone, not a new general lifting principle or a new matching
theorem.

At the documented check on 12 October 2026, the GitHub comparison returned
current `main` identical to the requested pin. The folder history showed only
its initial import; the complete upstream PR collection returned one already
merged PR, incorporated in the pin. No later relevant upstream change was
silently substituted. Bounded repository/web searches and the family-095
external discussion check did not locate an earlier smaller lift of this
fixed cone. **Not located is not verified originality.**

Full query provenance, retrieved evidence, pagination caveats, and the
distinction from the separate no-lift examples are recorded in
`audits/prior_work/PRIOR_WORK.md` and `evidence.json`. The pinned family-095
Lean scope concerns original-coordinate nonspectrahedrality of the degree-20
cone; this work does not claim a Lean formalization of either exact lift.

No publication, external communication, commit, push, or scheduled follow-up
was made. The result is a completed local research and verification handoff.

## Source references

Every path below is relative to the pinned folder
`preprints/An-Exact-Semidefinite-Lift-of-a-Nonspectrahedral-Hyperbolicity-Cone-October-5-2026/`.

- `build/main.tex`: title, date, stated dimensions, included sections.
- `build/sections/introduction.tex`: exact \(Q,\Phi,p,e,K\), 23-dimensional order,
  baseline theorem, relationship to the companion.
- `build/sections/angular.tex`: complex Gram factors and realification order.
- `build/sections/slices.tex`: skew blocks, \(\Phi\) identities, slice pencil,
  root identity and singular-slice equivalence.
- `build/sections/matching.tex`: unitary construction, both denominator branches,
  indefinite pairing, harmonic identification.
- `build/sections/lift.tex`: full functional construction, scalar certificate,
  range conditions, pseudoinverse completion, 100/307 coefficient rule.
- `build/references.bib`: cited mathematical background and attribution.

The unchanged files and exact source URLs are included in this archive, so the
proof and executable handoff can be read and replayed without network access.
