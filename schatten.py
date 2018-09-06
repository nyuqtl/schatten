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

def schatten(M, L, p) :
    dm = np.matrix(M, copy=True)
    for l in range(0, L) :
        dm[l, l] = 0.
    w, _ = np.linalg.eigh(dm)
    s = 0.
    for val in w :
        s += np.power(np.absolute(val), p)
    return np.power(s, 1./p)

def PLL(x, rho, L, n) :
    pn = np.zeros((L, L), dtype=np.complex128)
    for l1 in range(0, L) :
        for l2 in range(0, L) :
            if l1 != l2 :
                pn[l1, l2] = np.conj(x[l1, n])*x[l2, n]*rho[l1, l2]
    return pn

def instance(L, N) :
    x = np.zeros((L, N), dtype=np.complex128)
    # create X
    for l in range(0, L) :
        x[l,:] = np.matrix(qp.rand_ket(N).full()).T
    return x

def C2b(rho, x, L, N, p) :
    lhs = 0.
    for n in range(0, N) :
        pn = PLL(x, rho, L, n)
        lhs += schatten(pn, L, p)
    rhs = schatten(rho, L, p)
    return lhs <= rhs, lhs, rhs

dmhashes = []
check_duplicates = False

L = 4
N = 10
p = 1.
works = True
max = 1000000
for cnt in tqdm(range(0, max)) :
    seed = int(time.time()) + cnt
    np.random.seed(seed=seed)
    dm = np.matrix(qp.rand_dm(L).full())
    if check_duplicates :
        dmhashes.append(hash(str(dm)))
    x = instance(L, N)
    works, lhs, rhs = C2b(dm, x, L, N, p)
    if not works :
        print dm
        printerX(x, L, N)
        sys.exit(0)

if check_duplicates :
    print 'there were %s duplicating in density matrices' % str(len(dmhashes) - len(set(dmhashes)))

if works :
    print 'all %s random instances worked' % str(max)
