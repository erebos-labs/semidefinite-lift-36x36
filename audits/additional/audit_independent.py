if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

"""Third reconstruction. No imports from handoff programs or certificates."""
import sympy as s, json, pathlib, hashlib
ROOT=pathlib.Path(__file__).resolve().parents[2]
out=ROOT/'evidence/independent';out.mkdir(exist_ok=True)
w=s.Symbol('w',nonzero=True); I=s.I
c=(w+1/w)/2;sn=(w-1/w)/(2*I)
def star(v): return s.conjugate(v).subs(s.conjugate(w),1/w)
def re(v): return s.expand((v+star(v))/2)
def im(v): return s.expand((v-star(v))/(2*I))
U=s.Matrix([[-1/w,-w],[sn/w**2,sn*w**2],[c/w**2,-c*w**2]])/s.sqrt(2)
V=s.Matrix([[0,0],[I,I],[1,1]])/s.sqrt(2)
eta,xi=s.symbols('eta xi',real=True)
def B(r):
    a=s.zeros(4)
    for j in range(3):a[0,j+1]=r[j];a[j+1,0]=-r[j]
    return a
def sym(prefix):
    a=s.zeros(4)
    for i in range(4):
        for j in range(i,4):a[i,j]=a[j,i]=s.Symbol(prefix+str(i)+str(j),real=True)
    return a
X,Z=sym('x'),sym('z')
F=eta*U+xi*V
Rf=s.Matrix.hstack(*[v for j in range(2) for v in (F[:,j].applyfunc(re),F[:,j].applyfunc(im))])
L=s.diag(X,X,X,X,Z)
for j in range(4):
    bj=B(Rf[:,j]);L[4*j:4*j+4,16:20]=bj;L[16:20,4*j:4*j+4]=bj.T
P=s.zeros(20,36)
for group,(offset,phase) in enumerate([(0,1),(0,I),(0,w**-2),(0,I*w**-2),(8,1),(8,I),(8,w**2),(8,I*w**2)]):
    for j in range(4):P[offset+j,4*group+j]=re(s.sympify(phase));P[offset+4+j,4*group+j]=im(s.sympify(phase))
P[16:20,32:36]=s.eye(4)
poly=(P.T*L*P).applyfunc(s.expand)
pairs=[(i,j) for i in range(4) for j in range(i,4)]
coordinates=[]
for i,j in pairs:coordinates.extend([(X[i,j],k) for k in ('1','c2','s2')])
coordinates.extend((Z[i,j],'1') for i,j in pairs)
coordinates.extend([(eta,'c1'),(eta,'s1'),(eta,'c3'),(eta,'s3'),(xi,'1'),(xi,'c2'),(xi,'s2')])
variables=list(dict.fromkeys(a for a,k in coordinates))
lookup={key:k for k,key in enumerate(coordinates)}
mats=[s.zeros(36) for k in coordinates]
modechecks=0
for i in range(36):
    for j in range(36):
        for v in variables:
            e=s.expand(poly[i,j].coeff(v)); terms={}
            for t in s.Add.make_args(e):
                power=t.as_powers_dict().get(w,0)
                assert power==int(power)
                power=int(power);terms[power]=terms.get(power,0)+t/w**power
            modes={'1':terms.get(0,0)}
            for k in sorted({abs(t) for t in terms if t}):
                modes['c'+str(k)]=terms.get(k,0)+terms.get(-k,0)
                modes['s'+str(k)]=I*(terms.get(k,0)-terms.get(-k,0))
            for mode,a in modes.items():
                a=s.expand(a)
                if a!=0:
                    assert (v,mode) in lookup,(i,j,v,mode)
                    mats[lookup[v,mode]][i,j]=a
            modechecks+=1
payload=json.loads((ROOT/'coefficients.json').read_text())
for k,entries in enumerate(payload['coefficients']):
    expected=s.zeros(36)
    for i,j,(an,ad,bn,bd) in entries:expected[i,j]=expected[j,i]=s.Rational(an,ad)+s.Rational(bn,bd)*s.sqrt(2)
    assert mats[k]==expected,('full coefficient mismatch',k)
flat=s.Matrix([[m[i,j] for m in mats] for i,j in [(i,j) for i in range(36) for j in range(i,36)]])
outs=[lookup[X[i,j],'1'] for i,j in pairs]+[lookup[Z[i,j],'1'] for i,j in pairs]+[lookup[eta,'c1'],lookup[eta,'s1'],lookup[xi,'1']]
aux=[k for k in range(47) if k not in outs]
rankfull=flat.rank();rankaux=flat[:,aux].rank()
assert (rankfull,rankaux)==(47,24)
grams={name:sum((mats[k].T*mats[k] for k in indices),s.zeros(36)) for name,indices in [('full',range(47)),('output',outs),('auxiliary',aux)]}
kernels={name:g.nullspace() for name,g in grams.items()}
assert kernels['full']==[] and kernels['output']==[]
e33=s.eye(36)[:,33];assert kernels['auxiliary']==[e33]
auxrows=[pos for i,j in pairs for pos in ((i,8+j),(i,12+j))]+[(4,34),(0,34),(12,34),(12,35)]
minor=s.Matrix([[mats[k][i,j] for k in aux] for i,j in auxrows]); da=minor.det()
fullrows=[(i,j) for i,j in pairs]+[(32+i,32+j) for i,j in pairs]+auxrows+[(8,33),(12,33),(16,35)]
dm=s.Matrix([[mats[k][i,j] for k in range(47)] for i,j in fullrows]).det()
assert da==s.Rational(1,16) and dm==s.sqrt(2)/32
assert sum((mats[k] for k in [0,12,21,27,30,34,37,39]),s.zeros(36))==s.eye(36)
result={'status':'PASS','method':'direct symbolic P.T L P from angular factors; no handoff imports','sympy':s.__version__,'matrix_entries_compared':47*36*36,'directed_coordinate_contractions':modechecks,'nonzero_upper_entries':sum(sum(1 for i in range(36) for j in range(i,36) if m[i,j]!=0) for m in mats),'full_rank':rankfull,'aux_rank':rankaux,'common_kernels':{n:[list(v) for v in k] for n,k in kernels.items()},'aux_minor_det':str(da),'full_minor_det':str(dm)}
(out/'third-coefficients-ranks.json').write_text(json.dumps(result,default=str,indent=2)+'\n')
print(json.dumps(result,default=str,indent=2),flush=True)
