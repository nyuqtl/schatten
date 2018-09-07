import sys
import time
import qutip as qp
import numpy as np

from tqdm import tqdm

def printerX(x, L, N) :
    print '\nelements x that sum up to one over all n for each l'
    for l in range(0, L) :
        print ''
        total = 0.
        for n in range(0, N) :
            print 'l = %s, n = %s, x[%s, %s] = %s' % (str(l), str(n), str(l), str(n), str(x[l, n]))
            total += np.power(np.absolute(x[l, n]), 2.)
        print 'sums to %s' % str(total)

def schattenSingularSpec(M, L, p, rt=1.) :
    dm = np.zeros((L, L), dtype=np.float64)
    for l1 in range(0, L) :
        for l2 in range(0, L) :
            if l1 < l2 :
                # set upper right off-diagonal imaginary values to zero
                dm[l1, l2] = np.real(M[l1, l2])
            if l1 > l2 :
                # set lower left off-diagonal real values to zero
                dm[l1, l2] = np.imag(M[l1, l2])
    _, w, _ = np.linalg.svd(dm)
    s = 0.
    for val in w :
        s += np.power(np.absolute(val), p)
    return np.power(s, 1./(p*rt))

def schattenSingular(M, L, p, rt=1.) :
    dm = np.matrix(M, copy=True)
    for l in range(0, L) :
        dm[l, l] = 0.
    _, w, _ = np.linalg.svd(dm)
    s = 0.
    for val in w :
        s += np.power(np.absolute(val), p)
    return np.power(s, 1./(p*rt))

def schattenEigen(M, L, p, rt=1.) :
    dm = np.matrix(M, copy=True)
    for l in range(0, L) :
        dm[l, l] = 0.
    w, _ = np.linalg.eigh(dm)
    s = 0.
    for val in w :
        s += np.power(np.absolute(val), p)
    return np.power(s, 1./(p*rt))

def PLL(x, rho, L, n) :
    pn = np.zeros((L, L), dtype=np.complex128)
    for l1 in range(0, L) :
        for l2 in range(0, L) :
            if l1 != l2 :
                pn[l1, l2] = np.conj(x[l1, n])*x[l2, n]*rho[l1, l2]
    return np.matrix(pn)

def instance(L, N) :
    x = np.zeros((L, N), dtype=np.complex128)
    # create X
    for l in range(0, L) :
        x[l,:] = np.matrix(qp.rand_ket(N).full()).T
    return x

def C2b(rho, x, L, N, p, schatten) :
    lhs = 0.
    for n in range(0, N) :
        pn = PLL(x, rho, L, n)
        lhs += schatten(pn, L, p)
    rhs = schatten(rho, L, p)
    return lhs <= rhs, lhs, rhs

def C2bdagger(rho, x, L, N, p, schatten) :
    lhs = 0.
    for n in range(0, N) :
        pn = PLL(x, rho, L, n)
        pn = pn.getH()*pn
        lhs += schatten(pn, L, p, rt=2.)
    R = rho.getH()*rho
    rhs = schatten(R, L, p, rt=2.)
    return lhs <= rhs, lhs, rhs

def prnt(v) :
    if v :
        return str(1) + ' '
    return str(0) + ' '

def process(L, N, p, max, C2b, norm) :
    w = True
    for cnt in tqdm(range(0, max)) :
    #for cnt in range(0, max) :
        seed = int(time.time()) + cnt
        np.random.seed(seed=seed)
        dm = np.matrix(qp.rand_dm(L).full())
        x = instance(L, N)
        w, _, _ = C2b(dm, x, L, N, p, norm)
        #if False in [w1, w2, w3, w4, w5, w6] :
        if not w :
            break
    return w

dmhashes = []
check_duplicates = False

L = 2
N = 2
p = 2.
max = 100

w1 = process(L, N, p, max, C2b, schattenSingular)
w2 = process(L, N, p, max, C2b, schattenSingularSpec)
w3 = process(L, N, p, max, C2b, schattenEigen)
w4 = process(L, N, p, max, C2bdagger, schattenSingular)
w5 = process(L, N, p, max, C2bdagger, schattenSingularSpec)
w6 = process(L, N, p, max, C2bdagger, schattenEigen)

print '\nL = %s' % str(L)
print 'N = %s' % str(N)
print 'p = %s' % str(p)

print '\n+-- norm --+-- C2b --+-- C2b dagger --+'
print '| singular |   %s    |   %s           |' % (prnt(w1), prnt(w4))
print '| singspec |   %s    |   %s           |' % (prnt(w2), prnt(w5))
print '| eigenval |   %s    |   %s           |' % (prnt(w3), prnt(w6))
print '+----------+---------+----------------+\n'
