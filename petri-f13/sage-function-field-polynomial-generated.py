from sage.all import *
import json

Fp = GF(13)
Ru = PolynomialRing(Fp, "u")
u = Ru.gen()
g = u**6-u**5-3*u**4-3*u**3+u**2+5*u+4
k = GF(13**6, name="a", modulus=g)
a = k.gen()

def fv(value):
    coeffs = list(k(value).polynomial())
    coeffs += [Fp(0)] * (6-len(coeffs))
    return [int(coeffs[5-i]) for i in range(6)]

R4 = PolynomialRing(k, names=("x0","x1","x2","x3"))
x0,x1,x2,x3 = R4.gens()
Q = R4(((0)*x0*x0 + (0)*x0*x1 + (3*a**5 + 1*a**4 + 2*a**3 + 4*a**2 + 4*a + 2)*x0*x2 + (12)*x0*x3 + (10*a**5 + 12*a**4 + 11*a**3 + 9*a**2 + 9*a + 11)*x1*x1 + (1)*x1*x2 + (11*a**5 + 6*a**4 + 7*a**3 + 2*a**2 + 5*a + 10)*x1*x3 + (1)*x2*x2 + (3*a**5 + 8*a**4 + 3*a**3 + 4*a**2 + 2*a + 7)*x2*x3 + (7*a**5 + 3*a**4 + 3*a**3 + 8*a**2 + 11*a)*x3*x3))
C = R4(((0)*x0*x0*x0 + (0)*x0*x0*x1 + (4*a**5 + 10*a**4 + 2*a**3 + 5*a**2 + 8*a + 7)*x0*x0*x2 + (0)*x0*x0*x3 + (9*a**5 + 3*a**4 + 11*a**3 + 8*a**2 + 5*a + 6)*x0*x1*x1 + (6*a**5 + 2*a**4 + 1*a**3 + 10*a**2 + 5*a + 3)*x0*x1*x2 + (0)*x0*x1*x3 + (3*a**4 + 12*a**3 + 6*a**2 + 4*a + 2)*x0*x2*x2 + (0)*x0*x2*x3 + (0)*x0*x3*x3 + (7*a**5 + 11*a**4 + 12*a**3 + 3*a**2 + 8*a + 10)*x1*x1*x1 + (1)*x1*x1*x2 + (11*a**5 + 2*a**4 + 9*a**3 + 10*a**2 + 9*a + 4)*x1*x1*x3 + (10*a**4 + 12*a**3 + 6*a**2 + 6*a + 7)*x1*x2*x2 + (12*a**5 + 3*a**4 + 11*a**3 + 8*a**2 + 12*a + 8)*x1*x2*x3 + (5*a**5 + 3*a**3 + 3*a**2 + 7*a + 8)*x1*x3*x3 + (9*a**5 + 8*a**4 + 2*a**3 + 5*a**2 + 3*a + 7)*x2*x2*x2 + (7*a**5 + 4*a**4 + 11*a**3 + 1*a**2 + 10*a + 4)*x2*x2*x3 + (4*a**5 + 8*a**4 + 1*a**3 + 7*a**2 + 8*a + 12)*x2*x3*x3 + (11*a**5 + 4*a**4 + 10*a**3 + 6*a**2 + 9*a + 3)*x3*x3*x3))
N = R4(((1)*x0*x1*x3*x3*x3 + (1*a**5 + 5*a**4 + 7*a**3 + 7*a**2 + 5)*x0*x3*x3*x3*x3 + (4*a**5 + 5*a**4 + 5*a**3 + 7*a**2 + 3*a + 4)*x1*x1*x1*x2*x3 + (1*a**5 + 3*a**4 + 4*a**3 + 12*a + 5)*x1*x1*x1*x3*x3 + (9*a**5 + 8*a**4 + 8*a**3 + 6*a**2 + 10*a + 9)*x1*x1*x2*x2*x2 + (9*a**5 + 11*a**4 + 3*a**3 + 7*a**2 + 6*a + 7)*x1*x1*x2*x2*x3 + (2*a**4 + 10*a**3 + 3*a**2 + 4*a + 4)*x1*x1*x2*x3*x3 + (9*a**5 + 12*a**3 + 12*a**2 + 4*a + 6)*x1*x1*x3*x3*x3 + (6*a**5 + 2*a**4 + 9*a**3 + 8*a**2 + 1*a + 1)*x1*x2*x2*x2*x2 + (12*a**5 + 10*a**4 + 2*a**3 + 8*a**2 + 10*a + 6)*x1*x2*x2*x2*x3 + (12*a**5 + 6*a**4 + 7*a**3 + 3*a**2 + 4*a + 10)*x1*x2*x2*x3*x3 + (9*a**5 + 1*a**4 + 8*a**3 + 11*a**2 + 1)*x1*x2*x3*x3*x3 + (10*a**5 + 4*a**4 + 5*a**3 + 4*a**2 + 9*a + 4)*x1*x3*x3*x3*x3 + (7*a**5 + 5*a**4 + 10*a**3 + 11*a + 12)*x2*x2*x2*x2*x2 + (12*a**5 + 3*a**4 + 8*a**3 + 6*a**2 + 8*a)*x2*x2*x2*x2*x3 + (10*a**5 + 5*a**4 + 3*a**3 + 4*a**2 + 12*a + 1)*x2*x2*x2*x3*x3 + (7*a**5 + 5*a**4 + 12*a**2 + 5*a + 7)*x2*x2*x3*x3*x3 + (11*a**5 + 1*a**4 + 11*a**3 + 9*a**2 + 2*a + 4)*x2*x3*x3*x3*x3 + (1*a**5 + 2*a**4 + 7*a**3 + 6*a**2 + 12*a)*x3*x3*x3*x3*x3))
D = R4(((1)*x0*x0*x0*x0*x0 + (11*a**5 + 7*a**4 + 8*a**3 + 2*a**2 + 6*a + 4)*x0*x1*x1*x1*x1 + (1*a**5 + 5*a**4 + 7*a**3 + 5*a**2 + 10)*x0*x1*x3*x3*x3 + (10*a**5 + 1*a**4 + 3*a**3 + 2*a**2 + 5*a + 4)*x0*x3*x3*x3*x3 + (5*a**5 + 4*a**4 + 10*a**3 + 5*a**2 + 9*a)*x1*x1*x1*x1*x1 + (4*a**5 + 2*a**4 + 9*a**3 + 1*a**2 + 3*a + 6)*x1*x1*x1*x1*x2 + (7*a**5 + 3*a**4 + 11*a**3 + 8*a**2 + 7*a + 7)*x1*x1*x1*x2*x2 + (10*a**5 + 2*a**4 + 2*a**3 + 11*a**2 + 11*a + 6)*x1*x1*x1*x2*x3 + (8*a**5 + 3*a**4 + 3*a**2 + 9*a + 6)*x1*x1*x1*x3*x3 + (5*a**5 + 2*a**4 + 8*a**2 + 7*a + 8)*x1*x1*x2*x2*x2 + (3*a**5 + 11*a**4 + 6*a**3 + 1*a**2 + 12*a + 4)*x1*x1*x2*x2*x3 + (3*a**4 + 8*a**3 + 3*a**2 + 9*a)*x1*x1*x2*x3*x3 + (8*a**5 + 7*a**4 + 3*a**3 + 11*a**2 + 4*a)*x1*x1*x3*x3*x3 + (6*a**5 + 9*a**4 + 2*a**3 + 10*a**2 + 2*a + 8)*x1*x2*x2*x2*x2 + (11*a**5 + 2*a**4 + 12*a**3 + 12*a)*x1*x2*x2*x2*x3 + (7*a**5 + 3*a**4 + 5*a**3 + 12*a**2 + 4*a + 5)*x1*x2*x2*x3*x3 + (4*a**5 + 10*a**4 + 9*a**3 + 5*a**2 + 7*a + 12)*x1*x2*x3*x3*x3 + (12*a**5 + 9*a**4 + 12*a**3 + 10*a**2 + 3)*x1*x3*x3*x3*x3 + (4*a**5 + 3*a**4 + 5*a**3 + 11*a**2 + 7*a + 10)*x2*x2*x2*x2*x2 + (9*a**5 + 6*a**4 + 2*a**3 + 3*a**2 + 6*a + 12)*x2*x2*x2*x2*x3 + (5*a**5 + 6*a**4 + 8*a**3 + 10*a + 4)*x2*x2*x2*x3*x3 + (1*a**5 + 4*a**4 + 1*a**3 + 10*a**2 + 3)*x2*x2*x3*x3*x3 + (8*a**5 + 2*a**4 + 6*a**3 + 5*a**2 + 8*a + 10)*x2*x3*x3*x3*x3 + (4*a**5 + 5*a**4 + 8*a**3 + 11*a**2 + 10*a + 7)*x3*x3*x3*x3*x3))

Sx = PolynomialRing(k, names=("xx","y","z"))
xx,yx,zx = Sx.gens()
aff = R4.hom([xx,yx,zx,Sx(1)], Sx)
Qaff = aff(Q)
S = PolynomialRing(k, names=("y","z"))
y,z = S.gens()
drop = Sx.hom([S(0),y,z], S)
B = drop(Qaff.subs({xx:Sx(0)}))
A = drop(Qaff.subs({xx:Sx(1)}) - Qaff.subs({xx:Sx(0)}))
assert Qaff.degree(xx) == 1 and A != 0
FS = S.fraction_field()
sub = R4.hom([FS(-B/A),FS(y),FS(z),FS(1)], FS)

def substitute_clear(poly, power):
    value = sub(poly)
    return S(value * FS(A**power))

# The exact support shows degree_x0(C)<=2 and degree_x0(N),degree_x0(D)<=5.
assert C.degree(x0) <= 2 and N.degree(x0) <= 5 and D.degree(x0) <= 5
F = substitute_clear(C, 2)
Np = substitute_clear(N, 5)
Dp = substitute_clear(D, 5)
assert F != 0 and Np != 0 and Dp != 0

R = PolynomialRing(k, names=("tt","yy","zz"))
tt,yy,zz = R.gens()
FR = R(F(y=yy,z=zz))
NR = R(Np(y=yy,z=zz))
DR = R(Dp(y=yy,z=zz))
res = FR.resultant(NR-tt*DR, yy)
assert res != 0

Kt = FunctionField(k, "t")
t = Kt.gen()
P = PolynomialRing(Kt, "w")
w = P.gen()

def to_Kt(poly_t):
    # Coefficients are in k and tt is mapped to the function-field generator.
    return Kt(sum(k(coef) * t**exp[0] for exp,coef in poly_t.dict().items()))

resP = P([to_Kt(res.polynomial(zz)[i]) for i in range(res.degree(zz)+1)])
fac = resP.factor()
factors = []
for h,m in fac:
    hm = h.monic()
    factors.append({
        "degree": int(hm.degree()),
        "multiplicity": int(m),
        "coefficients": [
            {
                "numerator": [fv(v) for v in c.numerator().list()],
                "denominator": [fv(v) for v in c.denominator().list()],
            }
            for c in hm.list()
        ],
    })

out = {
    "status": "exact_function_field_elimination",
    "field": {"characteristic":13,"degree":6,"modulus":[1,12,10,10,1,5,4]},
    "chart": {
        "quadric_linear_coefficient_degree": [int(A.degree(y)),int(A.degree(z))],
        "quadric_remainder_degree": [int(B.degree(y)),int(B.degree(z))],
        "curve_degree": [int(F.degree(y)),int(F.degree(z))],
        "numerator_degree": [int(Np.degree(y)),int(Np.degree(z))],
        "denominator_degree": [int(Dp.degree(y)),int(Dp.degree(z))],
    },
    "raw_resultant": {
        "degree_t": int(res.degree(tt)),
        "degree_source": int(res.degree(zz)),
        "factor_count": len(factors),
    },
    "factors": factors,
}
out["degree23_factor_count"] = sum(1 for h in factors if h["degree"] == 23)
out["accepted"] = out["degree23_factor_count"] == 1
print("RESULT_JSON=" + json.dumps(out, sort_keys=True))
