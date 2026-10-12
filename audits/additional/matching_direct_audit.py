if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

"""Direct cleared-denominator generic-S matching, no handoff imports."""
import sympy as s, pathlib, json, time
started=time.monotonic();ROOT=pathlib.Path(__file__).resolve().parents[2]
w,v=s.symbols('w v',nonzero=True);I=s.I
def star(a):return s.conjugate(a).subs({s.conjugate(w):1/w,s.conjugate(v):1/v})
def real(a):return s.expand((a+star(a))/2)
def U(q):
    c=(q+1/q)/2;sn=(q-1/q)/(2*I)
    return s.Matrix([[-1/q,-q],[sn/q**2,sn*q**2],[c/q**2,-c*q**2]])/s.sqrt(2)
V=s.Matrix([[0,0],[I,I],[1,1]])/s.sqrt(2)
def Q(y):
    a,b,c=y;return s.Matrix([[a*a+b*b,-a*b,-a*c],[-a*b,b*b+c*c,-b*c],[-a*c,-b*c,c*c+a*a]])
def adj(a):return a.applyfunc(star).T
def check(a):
    for x in list(a) if isinstance(a,s.MatrixBase) else [a]:
        assert s.expand(x)==0
S=s.Matrix(3,3,lambda i,j:s.Symbol('S'+str(min(i,j))+str(max(i,j)),real=True))
tau,rho,eta,xi=s.symbols('tau rho eta xi',real=True)
R=tau*U(v)+rho*V;N=S*R;row=N[2,:]-I*N[1,:]
g=s.expand(row[0]+star(row[1]));h=s.expand(row[0]-star(row[1]))
delta=s.expand(g*star(g)+h*star(h));kappa=v**2/w**2
a=g*star(g)+h*star(h)*kappa;bb=g*h*(1-kappa)
Mn=s.Matrix([[a,star(bb)],[-bb,star(a)]]);C=(R*Mn).applyfunc(s.expand)
y=s.Matrix([eta*(w+1/w)/2,eta*(w-1/w)/(2*I),xi])
z=s.Matrix([tau*(v+1/v)/2,tau*(v-1/v)/(2*I),rho])
pairing=real(s.trace(adj(eta*U(w)+xi*V)*S*C))
polar=(Q(y+z)-Q(y)-Q(z))/2
check(pairing-delta*s.trace(S*polar))
def support(a):
    a=s.expand(a)
    return {int(term.as_powers_dict().get(w,0)) for term in s.Add.make_args(a) if term!=0}
supports=[sorted(set().union(*[support(C[i,j]) for i in range(3)])) for j in range(2)]
assert set(supports[0])<={0,-2} and set(supports[1])<={0,2}
# Verify scaled unitarity separately for fully generic complex g,h.
gr,gi,hr,hi=s.symbols('gr gi hr hi',real=True);gg=gr+I*gi;hh=hr+I*hi
dd=gg*star(gg)+hh*star(hh);aa=gg*star(gg)+hh*star(hh)*kappa;bbb=gg*hh*(1-kappa)
mm=s.Matrix([[aa,star(bbb)],[-bbb,star(aa)]])
check(mm*adj(mm)-dd**2*s.eye(2))
check(real(s.trace(adj(R)*S*R))-s.trace(S*Q(z)))
result={'status':'PASS','arithmetic':'exact cleared-denominator Laurent polynomials','generic_symmetric_S_no_definiteness_assumption':True,'direct_matching_pairing_numerator_zero':True,'C_numerator_column_supports':supports,'generic_rotation_scaled_unitary':True,'delta_zero_branch':'Reviewed separately in written proof and reproduced 42-check replay; no division at zero','seconds_for_information':time.monotonic()-started}
(ROOT/'evidence/independent/direct-generic-matching.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
