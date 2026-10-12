if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

import sympy as s, pathlib, json
ROOT=pathlib.Path(__file__).resolve().parents[2]; out=ROOT/'evidence/independent'
t=s.Symbol('t');rt=s.sqrt(2)
def Q(y):
    a,b,c=y;return s.Matrix([[a*a+b*b,-a*b,-a*c],[-a*b,b*b+c*c,-b*c],[-a*c,-b*c,c*c+a*a]])
def phi(W,y):
    q=Q(y);return s.Matrix.vstack(s.Matrix.hstack(s.Matrix([[s.trace(q*W[1:,1:])]]),-W[0:1,1:]*q),s.Matrix.hstack(-q*W[1:,0:1],W[0,0]*q))
def B(r):return s.Matrix([[0,*r],[-r[0],0,0,0],[-r[1],0,0,0],[-r[2],0,0,0]])
def realify(F):return s.Matrix.hstack(*[v for j in range(2) for v in (s.re(F[:,j]),s.im(F[:,j]))])
V=s.Matrix([[0,0],[s.I,s.I],[1,1]])/rt
bs=[B(realify(V)[:,j]) for j in range(4)]
def slice(X,Z):
    b=s.Matrix.vstack(*bs);return s.Matrix.vstack(s.Matrix.hstack(s.diag(X,X,X,X),b),s.Matrix.hstack(b.T,Z))
def root_identity(X,Z,y,L,label):
    Xt=t*s.eye(4)-X;Zt=t*s.eye(4)-Z
    original=(Xt.det()*Zt-phi(Xt.adjugate(),-y)).det(method='domain-ge')
    cp=L.charpoly(t).as_expr();assert s.expand(original-cp)==0
    poly=s.Poly(cp,t);assert poly.degree()==20 and poly.LC()==1
    multiplicity=next(k for k in range(21) if poly.nth(k)!=0)
    return {'case':label,'degree':20,'original_p_equals_characteristic_polynomial':True,'zero_root_multiplicity':multiplicity,'factorization':str(s.factor(cp))}
results=[]
y=s.Matrix([0,0,1]);X=s.Matrix([[2,0,1,0],[0,0,0,0],[1,0,1,0],[0,0,0,3]])
# Construct the Moore-Penrose inverse from the active block, not a numeric pinv.
idx=[0,2,3];Xp=s.zeros(4);active=X.extract(idx,idx).inv()
for i,a in enumerate(idx):
    for j,b in enumerate(idx):Xp[a,b]=active[i,j]
assert Xp*X*Xp==Xp and X*Xp==s.diag(1,0,1,1)
for bj in bs:assert X*Xp*bj==bj
residual=s.diag(0,1,0,0);Z=phi(Xp,y)+residual;L=slice(X,Z)
C=s.Matrix.hstack(s.eye(16),s.Matrix.vstack(*[Xp*bj for bj in bs]))
assert L==C.T*s.diag(X,X,X,X)*C+s.diag(s.zeros(16),residual)
assert X.extract(idx,idx)[:1,:1].det()>0 and X.extract(idx,idx)[:2,:2].det()>0 and X.extract(idx,idx).det()>0
assert L.rank()==13
results.append(root_identity(X,Z,y,L,'off-diagonal-singular-X-with-positive-residual'))
# Non-cardinal angle verifies every sin2 placement in the finite witness.
w=s.Rational(3,5)+s.Rational(4,5)*s.I;c2=s.re(w**2);s2=s.im(w**2)
P=s.zeros(20,36)
for group,(offset,phase) in enumerate([(0,1),(0,s.I),(0,w**-2),(0,s.I*w**-2),(8,1),(8,s.I),(8,w**2),(8,s.I*w**2)]):
    for j in range(4):P[offset+j,4*group+j]=s.re(phase);P[offset+4+j,4*group+j]=s.im(phase)
P[16:20,32:36]=s.eye(4)
pairs=[(i,j) for i in range(4) for j in range(i,4)]
moments=[v for i,j in pairs for v in (X[i,j],X[i,j]*c2,X[i,j]*s2)]+[Z[i,j] for i,j in pairs]+[0,0,0,0,1,c2,s2]
payload=json.loads((ROOT/'coefficients.json').read_text());H=s.zeros(36)
for k,entries in enumerate(payload['coefficients']):
    for i,j,(an,ad,bn,bd) in entries:
        a=moments[k]*(s.Rational(an,ad)+s.Rational(bn,bd)*rt);H[i,j]+=a
        if i!=j:H[j,i]+=a
assert (H-P.T*L*P).applyfunc(s.simplify)==s.zeros(36)
assert H.rank()==13
results.append({'case':'noncardinal-angle-finite-witness','w':'3/5+4i/5','cos2':str(c2),'sin2':str(s2),'slice_rank':13,'lift_rank':13,'H_equals_PtLP':True})
Xbad=s.diag(1,0,0,0);Zbad=phi(Xbad,y);Lbad=slice(Xbad,Zbad)
assert Zbad==s.diag(0,0,1,1)
assert Lbad.extract([3,16],[3,16]).det()==-s.Rational(1,2)
assert Xbad*bs[0]!=bs[0]
results.append(root_identity(Xbad,Zbad,y,Lbad,'pseudoinverse-only-invalid-point'))
results.append({'case':'pseudoinverse-only-invalid-point','Schur_residual':'zero','negative_principal_minor':'-1/2','range_condition_fails':True})
# At X=0 the Schur residual test alone says nothing about nonzero y.
Xzero=s.zeros(4);Zzero=s.eye(4);Lzero=slice(Xzero,Zzero)
assert Lzero.extract([3,16],[3,16]).det()==-s.Rational(1,2)
results.append(root_identity(Xzero,Zzero,y,Lzero,'X-zero-nonzero-y-invalid-point'))
# Check multiplicities even where the polynomial vanishes at the point itself.
results.append(root_identity(s.zeros(4),s.zeros(4),s.zeros(3,1),s.zeros(20),'origin'))
result={'status':'PASS','arithmetic':'exact symbolic/rational','checks':results,'limits':'Finite regressions supplement the universal written singular proof.'}
(out/'singular-root-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
