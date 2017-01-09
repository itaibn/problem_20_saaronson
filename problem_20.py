"""
Creates distributions D0,D1 where you can sampe (D0^2 + D1^2) easily but neither
D0 nor D1.

https://itaibn.wordpress.com/2013/09/29/a-solution-to-problem-20-in-scott-aaronsons-projects-aplenty/

NOTE: This solution was later found to be invalid, see blog post for details.
"""

from sage.all import *

def symbol(x, R, r):
    ring = IntegerModRing(R)
    return ring(x) ** ((R-1) / r)

def binomial_distribution(n, p):
    """
    Sample from binomial distribution B(n, p). Currently samples a corresponding
    Gaussian, which is close enough for the purpose of this program.
    """
    mean = n * p
    variance = n * p * (1-p)
    D = RealDistribution('gaussian', sqrt(variance))
    return floor(D.get_random_element() + mean)

def gen_A(l):
    """
    Generate DistributionData instance (see below)
    """
    size = 2 ** l
    p = random_prime(size)
    q0 = random_prime(size)
    q1 = random_prime(size)

    def random_prime_1_mod_m(limit, m):
        while True:
            p = randint(1, limit)*m + 1
            if is_prime(p):
                return p

    P = random_prime_1_mod_m(size, p)
    Q = random_prime_1_mod_m(size, q0*q1)
    N = P*Q

    ZP = IntegerModRing(P)
    ZQ = IntegerModRing(Q)
    ZN = IntegerModRing(N)

    n = l*4
    m = size
    gp = ZP.random_element()**2
    g = symbol(gp, P, p)

    # A canonical non-square with Jacobi symbol 1. Is also a p*q0*q1'th root for
    # convenience.
    nsP = ZP(primitive_root(P)) ** p
    nsQ = ZQ(primitive_root(Q)) ** (q0 * q1)
    ns = ZN(CRT(ZZ(nsP), ZZ(nsQ), P, Q))

    u, up, r, rp = [], [], [], []
    for _ in range(n):
        up.append(ZQ.random_element() ** (2*q1))
        rp.append(ZQ.random_element() ** (2*q0))
        u.append(symbol(up[-1], Q, q0))
        r.append(symbol(up[-1], Q, q1))

    A = matrix(ZN, n)
    b = matrix(n) # bit-matrix for weather something is a square

    for i in range(n):
        for j in range(n):
            # Generate A[i,j]
            if i > j:
                s = b[j, i]
            else:
                s = randint(0, 1)
            b[i, j] = s
            # Generate the value mod P and mod Q seperately as aP and aQ
            # respectively.
            gexp = binomial_distribution(m, 0.5)
            aP = gp ** gexp * ZP.random_element() ** (2*p)
            aQ = rp[i] * up[j] * ZQ.random_element() ** (2*q0*q1)
            A[i, j] = CRT(ZZ(aP), ZZ(aQ), P, Q) * ns ** s

    return DistributionData(N, n, A)

class DistributionData():
    """
    The data needed to sample (D0^2+D1^2)/2.
    """
    def __init__(self, N, n, A):
        self.N = N
        self.n = n
        self.A = A

    def gen_pair(self):
        Sn = Permutations(self.n)
        pi = Sn.random_element()
        return self.f(pi), self.f(pi.inverse())

    def f(self, pi):
        entries = [self.A[i, pi(i+1)-1] for i in range(self.n)]
        return prod(entries)
