from sage import *

def symbol(x, R, r):
    ring = IntegerModRing(R)
    return ring(x) ^ ((R-1) / r)

def gen_A(l):
    size = 2^l
    p = random_prime(size)
    q0 = random_prime(size)
    q1 = random_prime(size)

    def random_prime_1_mod_m(limit, m):
        while True:
            p = random(1, limit)*m + 1
            if is_prime(p):
                return p

    P = random_prime_1_mod_m(size, p)
    Q = random_prime_1_mod_m(size, q0*q1)
    N = P*Q

    ZP = IntegerModRing(P)
    ZQ = IntegerModRing(Q)
    ZN = IntegerModRing(N)

    n = l*4
    gp = ZP.random_element()
    g = symbol(gp, P, p)
    u, up, r, rp = []
    for _ in range(n):
        up.append(ZQ.random_element()^q1)
        rp.append(ZQ.random_element()^q0)
        u.append(symbol(up[-1], Q, q0))
        r.append(symbol(up[-1], Q, q1))
    A = matrix(n)

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
