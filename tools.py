import numpy as np
import qutip as qp

def testDM(rho, silent=True) :
    w, _ = np.linalg.eigh(rho)
    trace = np.trace(rho)
    t1 = np.isclose([trace], [1.])[0]
    t2 = (w >= 0.).all()
    if not t2 :
        print w
    t3 = (rho.H == rho).all()
    if not silent :
        print 'trace is one: %s' % str(t1)
        print 'eigenvalues positive: %s' % str(t2)
        print 'is Hermitian: %s' % str(t3)
    return t1 and t2 and t3

def schattenP(M, L, p, rt=1.) :
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

def C2b(rho, x, L, N, p) :
    lhs = 0.
    for n in range(0, N) :
        pn = PLL(x, rho, L, n)
        lhs += schattenP(pn, L, p)
    rhs = schattenP(rho, L, p)
    return lhs, rhs

def instanceX(L, N) :
    x = np.zeros((L, N), dtype=np.complex128)
    # create X
    for l in range(0, L) :
        x[l,:] = np.matrix(qp.rand_ket(N).full()).T
    return x

def Bn(M, L) :
    dm = np.zeros((L, L), dtype=np.float64)
    for l1 in range(0, L) :
        for l2 in range(0, L) :
            if l1 < l2 :
                # set upper right off-diagonal imaginary values to zero
                dm[l1, l2] = np.real(M[l1, l2])
            if l1 > l2 :
                # set lower left off-diagonal real values to zero
                dm[l1, l2] = np.imag(M[l1, l2])
    return np.matrix(dm)

def BiB(rho, x, L, N, p) :
    lhs = 0.
    for n in range(0, N) :
        pn = Bn(PLL(x, rho, L, n), L)
        pn = pn + 1j*pn.getH()
        lhs += schattenP(pn, L, p)
    B = Bn(rho, L)
    R = B+1j*B.getH()
    rhs = schattenP(R, L, p)
    return lhs, rhs
