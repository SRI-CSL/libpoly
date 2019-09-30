import polypy
import random
import sympy

PASS = 0
FAIL = 0

"""
 Call to start a unit test.
"""
def start(description):
    print("\t* Checking {0}".format(description))

"""
 Call to check a unit test.
"""
def check(ok):
    global PASS, FAIL
    if ok == True:
        PASS = PASS + 1
    elif ok == False:
        FAIL = FAIL + 1
    else:
        print("Didn't get True/False")

"""
 Make a random polynomial (degree bound and M bound on the coefficient
 magnitude.
"""
def random_upolynomial(K, degree, M, lc = None):
    coeff = [random.randint(-M, M) for _ in range(degree)]
    if (lc is None):
        coeff.append(random.randint(1, M))
    else:
        coeff.append(lc)
    p = polypy.UPolynomial(K, coeff)
    return p

def random_monomial(degree, M, p_vars):
    m = random.randint(1, M)
    if random.randint(0, 1):
        m_vars = [p_vars[i] for i in sorted(random.sample(range(len(p_vars)), random.randint(0, len(p_vars))))]
        for var in m_vars:
            deg = random.randint(1, degree)
            m = m*(var**deg)
    return m

"""
 Make a random polynomial (degree bound and M bound on the coefficients).
"""
def random_polynomial(degree, M, p_vars, trials):
    # Generate monomials
    p = random.randint(1, M)
    for _ in range(trials):
        p = p + random_monomial(degree, M, p_vars)
    while (not isinstance(p, polypy.Polynomial)) or (p.degree() == 0):
        p = p + random_monomial(degree, M, p_vars)
    return p

class SympyWrapper:

    enabled = True

    def sympy_from_upolynomial(self, p):
        coeffs = p.coefficients()
        sympy_p = 0
        x = sympy.symbols('x')
        for d, c in enumerate(coeffs):
            sympy_p = sympy_p + c*x**d
        sympy_p = sympy.simplify(sympy_p)
        return sympy_p

    def sympy_factor(self, p):
        sympy_p = self.sympy_from_upolynomial(p)
        if (p.ring().modulus() is None):
            return sympy.factor_list(sympy_p)
        else:
            return sympy.factor_list(sympy_p, modulus=p.ring().modulus())

    def sympy_gcd(self, p, q):
        sympy_p = self.sympy_from_upolynomial(p)
        sympy_q = self.sympy_from_upolynomial(q)
        if (p.ring().modulus() is None):
            return sympy.gcd(sympy_p, sympy_q)
        else:
            return sympy.gcd(sympy_p, sympy_q, modulus=p.ring().modulus())

    def sympy_extended_gcd(self, p, q):
        sympy_p = self.sympy_from_upolynomial(p)
        sympy_q = self.sympy_from_upolynomial(q)
        if (p.ring().modulus() is None):
            (u, v, d) = sympy.gcdex(sympy_p, sympy_q)
        else:
            (u, v, d) = sympy.gcdex(sympy_p, sympy_q, modulus=p.ring().modulus())
        return (d.simplify(), u.simplify(), v.simplify())

    def sympy_factor_count(self, p):
        factors = self.sympy_factor(p)
        return len(factors[1])

    def check_gcd(self, p, q, gcd):
        if not self.enabled:
            return True
        gcd_gold = self.sympy_gcd(p, q).simplify()
        ok = self.sympy_from_upolynomial(gcd) == gcd_gold
        if (not ok):
            print("Wrong gcd")
            print("p = {0}".format(p))
            print("q = {0}".format(q))
            print("gcd = {0}".format(gcd))
            print("expected = {0}".format(gcd_gold))
        return ok

    def check_extended_gcd(self, p, q, gcd, u, v):
        if not self.enabled:
            return True;
        (sympy_gcd, sympy_u, sympy_v) = self.sympy_extended_gcd(p, q)
        ok = self.sympy_from_upolynomial(gcd) == sympy_gcd
        if (not ok):
            print("Wrong gcd")
            print("p = {0}".format(p))
            print("q = {0}".format(q))
            print("gcd = {0}".format(gcd))
            print("expected = {0}".format(sympy_gcd))
            return False
        ok = self.sympy_from_upolynomial(u) == sympy_u
        if (not ok):
            print("Wrong u")
            print("p = {0}".format(p))
            print("q = {0}".format(q))
            print("u = {0}".format(u))
            print("expected = {0}".format(sympy_u))
            return False
        ok = self.sympy_from_upolynomial(v) == sympy_v
        if (not ok):
            print("Wrong v")
            print("p = {0}".format(p))
            print("q = {0}".format(q))
            print("v = {0}".format(v))
            print("expected = {0}".format(sympy_v))
            return False
        return True


"""
 By default sympy is enabled.
"""
sympy_checker = SympyWrapper();

"""
 Initialize the testing.
"""
def init():
    global PASS, FAIL
    PASS = 0
    FAIL = 0
    random.seed(0)

import sys

"""
 Disable output buffering.
"""
class NoBufferWrapper(object):
    def __init__(self, stream):
        self.stream = stream
    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)

sys.stdout = NoBufferWrapper(sys.stdout)
