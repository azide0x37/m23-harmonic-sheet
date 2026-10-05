from sage.all import *
import sage.version
import json
import time

started = time.monotonic()
Fp = GF(13)
Rz = PolynomialRing(Fp, "z")
z = Rz.gen()
g = z**6-z**5-3*z**4-3*z**3+z**2+5*z+4
k = GF(13**6, name="a", modulus=g)
a = k.gen()

def field_vector(value):
    coeffs = list(value.polynomial())
    coeffs += [Fp(0)] * (6-len(coeffs))
    return [int(coeffs[5-i]) for i in range(6)]

def hp_text(ideal):
    return str(ideal.hilbert_polynomial())

R4 = PolynomialRing(k, names=("x0","x1","x2","x3"), order="degrevlex")
x0,x1,x2,x3 = R4.gens()
Q = R4(((0)*x0*x0 + (0)*x0*x1 + (3*a**5 + 1*a**4 + 2*a**3 + 4*a**2 + 4*a + 2)*x0*x2 + (12)*x0*x3 + (10*a**5 + 12*a**4 + 11*a**3 + 9*a**2 + 9*a + 11)*x1*x1 + (1)*x1*x2 + (11*a**5 + 6*a**4 + 7*a**3 + 2*a**2 + 5*a + 10)*x1*x3 + (1)*x2*x2 + (3*a**5 + 8*a**4 + 3*a**3 + 4*a**2 + 2*a + 7)*x2*x3 + (7*a**5 + 3*a**4 + 3*a**3 + 8*a**2 + 11*a)*x3*x3))
C = R4(((0)*x0*x0*x0 + (0)*x0*x0*x1 + (4*a**5 + 10*a**4 + 2*a**3 + 5*a**2 + 8*a + 7)*x0*x0*x2 + (0)*x0*x0*x3 + (9*a**5 + 3*a**4 + 11*a**3 + 8*a**2 + 5*a + 6)*x0*x1*x1 + (6*a**5 + 2*a**4 + 1*a**3 + 10*a**2 + 5*a + 3)*x0*x1*x2 + (0)*x0*x1*x3 + (3*a**4 + 12*a**3 + 6*a**2 + 4*a + 2)*x0*x2*x2 + (0)*x0*x2*x3 + (0)*x0*x3*x3 + (7*a**5 + 11*a**4 + 12*a**3 + 3*a**2 + 8*a + 10)*x1*x1*x1 + (1)*x1*x1*x2 + (11*a**5 + 2*a**4 + 9*a**3 + 10*a**2 + 9*a + 4)*x1*x1*x3 + (10*a**4 + 12*a**3 + 6*a**2 + 6*a + 7)*x1*x2*x2 + (12*a**5 + 3*a**4 + 11*a**3 + 8*a**2 + 12*a + 8)*x1*x2*x3 + (5*a**5 + 3*a**3 + 3*a**2 + 7*a + 8)*x1*x3*x3 + (9*a**5 + 8*a**4 + 2*a**3 + 5*a**2 + 3*a + 7)*x2*x2*x2 + (7*a**5 + 4*a**4 + 11*a**3 + 1*a**2 + 10*a + 4)*x2*x2*x3 + (4*a**5 + 8*a**4 + 1*a**3 + 7*a**2 + 8*a + 12)*x2*x3*x3 + (11*a**5 + 4*a**4 + 10*a**3 + 6*a**2 + 9*a + 3)*x3*x3*x3))
N = R4(((1)*x0*x1*x3*x3*x3 + (1*a**5 + 5*a**4 + 7*a**3 + 7*a**2 + 5)*x0*x3*x3*x3*x3 + (4*a**5 + 5*a**4 + 5*a**3 + 7*a**2 + 3*a + 4)*x1*x1*x1*x2*x3 + (1*a**5 + 3*a**4 + 4*a**3 + 12*a + 5)*x1*x1*x1*x3*x3 + (9*a**5 + 8*a**4 + 8*a**3 + 6*a**2 + 10*a + 9)*x1*x1*x2*x2*x2 + (9*a**5 + 11*a**4 + 3*a**3 + 7*a**2 + 6*a + 7)*x1*x1*x2*x2*x3 + (2*a**4 + 10*a**3 + 3*a**2 + 4*a + 4)*x1*x1*x2*x3*x3 + (9*a**5 + 12*a**3 + 12*a**2 + 4*a + 6)*x1*x1*x3*x3*x3 + (6*a**5 + 2*a**4 + 9*a**3 + 8*a**2 + 1*a + 1)*x1*x2*x2*x2*x2 + (12*a**5 + 10*a**4 + 2*a**3 + 8*a**2 + 10*a + 6)*x1*x2*x2*x2*x3 + (12*a**5 + 6*a**4 + 7*a**3 + 3*a**2 + 4*a + 10)*x1*x2*x2*x3*x3 + (9*a**5 + 1*a**4 + 8*a**3 + 11*a**2 + 1)*x1*x2*x3*x3*x3 + (10*a**5 + 4*a**4 + 5*a**3 + 4*a**2 + 9*a + 4)*x1*x3*x3*x3*x3 + (7*a**5 + 5*a**4 + 10*a**3 + 11*a + 12)*x2*x2*x2*x2*x2 + (12*a**5 + 3*a**4 + 8*a**3 + 6*a**2 + 8*a)*x2*x2*x2*x2*x3 + (10*a**5 + 5*a**4 + 3*a**3 + 4*a**2 + 12*a + 1)*x2*x2*x2*x3*x3 + (7*a**5 + 5*a**4 + 12*a**2 + 5*a + 7)*x2*x2*x3*x3*x3 + (11*a**5 + 1*a**4 + 11*a**3 + 9*a**2 + 2*a + 4)*x2*x3*x3*x3*x3 + (1*a**5 + 2*a**4 + 7*a**3 + 6*a**2 + 12*a)*x3*x3*x3*x3*x3))
D = R4(((1)*x0*x0*x0*x0*x0 + (11*a**5 + 7*a**4 + 8*a**3 + 2*a**2 + 6*a + 4)*x0*x1*x1*x1*x1 + (1*a**5 + 5*a**4 + 7*a**3 + 5*a**2 + 10)*x0*x1*x3*x3*x3 + (10*a**5 + 1*a**4 + 3*a**3 + 2*a**2 + 5*a + 4)*x0*x3*x3*x3*x3 + (5*a**5 + 4*a**4 + 10*a**3 + 5*a**2 + 9*a)*x1*x1*x1*x1*x1 + (4*a**5 + 2*a**4 + 9*a**3 + 1*a**2 + 3*a + 6)*x1*x1*x1*x1*x2 + (7*a**5 + 3*a**4 + 11*a**3 + 8*a**2 + 7*a + 7)*x1*x1*x1*x2*x2 + (10*a**5 + 2*a**4 + 2*a**3 + 11*a**2 + 11*a + 6)*x1*x1*x1*x2*x3 + (8*a**5 + 3*a**4 + 3*a**2 + 9*a + 6)*x1*x1*x1*x3*x3 + (5*a**5 + 2*a**4 + 8*a**2 + 7*a + 8)*x1*x1*x2*x2*x2 + (3*a**5 + 11*a**4 + 6*a**3 + 1*a**2 + 12*a + 4)*x1*x1*x2*x2*x3 + (3*a**4 + 8*a**3 + 3*a**2 + 9*a)*x1*x1*x2*x3*x3 + (8*a**5 + 7*a**4 + 3*a**3 + 11*a**2 + 4*a)*x1*x1*x3*x3*x3 + (6*a**5 + 9*a**4 + 2*a**3 + 10*a**2 + 2*a + 8)*x1*x2*x2*x2*x2 + (11*a**5 + 2*a**4 + 12*a**3 + 12*a)*x1*x2*x2*x2*x3 + (7*a**5 + 3*a**4 + 5*a**3 + 12*a**2 + 4*a + 5)*x1*x2*x2*x3*x3 + (4*a**5 + 10*a**4 + 9*a**3 + 5*a**2 + 7*a + 12)*x1*x2*x3*x3*x3 + (12*a**5 + 9*a**4 + 12*a**3 + 10*a**2 + 3)*x1*x3*x3*x3*x3 + (4*a**5 + 3*a**4 + 5*a**3 + 11*a**2 + 7*a + 10)*x2*x2*x2*x2*x2 + (9*a**5 + 6*a**4 + 2*a**3 + 3*a**2 + 6*a + 12)*x2*x2*x2*x2*x3 + (5*a**5 + 6*a**4 + 8*a**3 + 10*a + 4)*x2*x2*x2*x3*x3 + (1*a**5 + 4*a**4 + 1*a**3 + 10*a**2 + 3)*x2*x2*x3*x3*x3 + (8*a**5 + 2*a**4 + 6*a**3 + 5*a**2 + 8*a + 10)*x2*x3*x3*x3*x3 + (4*a**5 + 5*a**4 + 8*a**3 + 11*a**2 + 10*a + 7)*x3*x3*x3*x3*x3))
I = R4.ideal(Q,C)
irrelevant = R4.ideal(x0,x1,x2,x3)
print("curve ideal ready", flush=True)

Braw = I + R4.ideal(N,D)
B, Bexp = Braw.saturation(irrelevant)
print("base saturation", Bexp, hp_text(B), flush=True)

JN, _ = (I + R4.ideal(N)).saturation(irrelevant)
JD, _ = (I + R4.ideal(D)).saturation(irrelevant)
ZN, ZNexp = JN.saturation(B)
ZD, ZDexp = JD.saturation(B)
ZN, _ = ZN.saturation(irrelevant)
ZD, _ = ZD.saturation(irrelevant)
print("zero/pole residual", hp_text(ZN), hp_text(ZD), flush=True)

# Ramification at finite nonzero target values.  The 3x3 minors say that the
# pencil section N-tD is tangent to the smooth complete intersection.
R5 = PolynomialRing(k, names=("x0","x1","x2","x3","t"), order="degrevlex")
x0,x1,x2,x3,t = R5.gens()
Q5,C5,N5,D5 = map(R5,(Q,C,N,D))
F5 = N5-t*D5
variables = (x0,x1,x2,x3)
rows = ([Q5.derivative(v) for v in variables],
        [C5.derivative(v) for v in variables],
        [F5.derivative(v) for v in variables])
minors = []
for omitted in range(4):
    cols = [j for j in range(4) if j != omitted]
    minors.append(matrix(R5, [[rows[i][j] for j in cols] for i in range(3)]).det())
J = R5.ideal(Q5,C5,F5,*minors)
m5 = R5.ideal(x0,x1,x2,x3)
B5 = R5.ideal(Q5,C5,N5,D5)
J, Jirr_exp = J.saturation(m5)
print("ramification irrelevant saturation", Jirr_exp, flush=True)
J, Jbase_exp = J.saturation(B5)
print("ramification base saturation", Jbase_exp, flush=True)
J, Jzero_exp = J.saturation(R5.ideal(t))
print("ramification zero-value saturation", Jzero_exp, flush=True)
E = J.elimination_ideal([x0,x1,x2,x3])
print("elimination", E, flush=True)
egens = [p for p in E.gens() if p != 0]

result = {
  "status": "sage_exact_passport_probe",
  "sage_version": sage.version.version,
  "field": {"characteristic": 13, "degree": 6},
  "curve_ideal_dimension": int(I.dimension()),
  "base": {
    "saturation_exponent": int(Bexp),
    "dimension": int(B.dimension()),
    "hilbert_polynomial": hp_text(B),
    "groebner_basis_size": len(B.groebner_basis()),
  },
  "zero_fibre_after_base_removal": {
    "saturation_exponent": int(ZNexp),
    "dimension": int(ZN.dimension()),
    "hilbert_polynomial": hp_text(ZN),
    "radical_hilbert_polynomial": hp_text(ZN.radical()),
  },
  "pole_fibre_after_base_removal": {
    "saturation_exponent": int(ZDexp),
    "dimension": int(ZD.dimension()),
    "hilbert_polynomial": hp_text(ZD),
    "radical_hilbert_polynomial": hp_text(ZD.radical()),
  },
  "finite_nonzero_ramification": {
    "irrelevant_saturation_exponent": int(Jirr_exp),
    "base_saturation_exponent": int(Jbase_exp),
    "zero_value_saturation_exponent": int(Jzero_exp),
    "ideal_dimension": int(J.dimension()),
    "elimination_generators": [str(p) for p in egens],
  },
  "elapsed_seconds": time.monotonic()-started,
}
if len(egens) == 1 and egens[0].degree(t) == 1:
    poly = egens[0]
    lam = k(-poly.subs({t:0}) / poly.monomial_coefficient(t))
    result["finite_nonzero_ramification"]["branch_value"] = field_vector(lam)
    result["finite_nonzero_ramification"]["branch_value_text"] = str(lam)
    Slambda = N-lam*D
    Flambda, _ = (I + R4.ideal(Slambda)).saturation(irrelevant)
    Fibre, fibre_base_exp = Flambda.saturation(B)
    Fibre, _ = Fibre.saturation(irrelevant)
    Radical = Fibre.radical()
    radical_square = I + R4.ideal([
        left*right for left in Radical.gens() for right in Radical.gens()
    ])
    square_contained = bool(radical_square <= Fibre)
    fibre_hp = Fibre.hilbert_polynomial()
    radical_hp = Radical.hilbert_polynomial()
    result["third_fibre"] = {
      "base_saturation_exponent": int(fibre_base_exp),
      "dimension": int(Fibre.dimension()),
      "hilbert_polynomial": str(fibre_hp),
      "radical_hilbert_polynomial": str(radical_hp),
      "radical_square_contained_in_fibre_ideal": square_contained,
      "fibre_ideal_contained_in_radical": bool(Fibre <= Radical),
      "length": int(fibre_hp),
      "radical_length": int(radical_hp),
      "excess_length": int(fibre_hp-radical_hp),
    }
    if int(fibre_hp) == 23 and int(radical_hp) == 15 and square_contained:
        result["third_fibre"].update({
          "multiplicity_partition_geometric": [2]*8 + [1]*7,
          "ramification_contribution": 8,
          "passport_certified": True,
        })
print("RESULT_JSON=" + json.dumps(result, sort_keys=True), flush=True)
