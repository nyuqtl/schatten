import numpy as np
import qutip as qp

def instanceX(L, N) :
    x = np.zeros((L, N), dtype=np.complex128)
    # create X
    for l in range(0, L) :
        x[l,:] = np.matrix(qp.rand_ket(N).full()).T
    return x

def PLL(x, rho, L, n) :
    pn = np.zeros((L, L), dtype=np.complex128)
    for l1 in range(0, L) :
        for l2 in range(0, L) :
            if l1 != l2 :
                pn[l1, l2] = np.conj(x[l1, n])*x[l2, n]*rho[l1, l2]
    return np.matrix(pn)

def wrapperBiB(rho, L) :
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
    B = Bn(rho, L)
    return np.exp(-1j*np.pi/4.)*B + np.exp(1j*np.pi/4.)*B.getH()

def wrapperIdentity(rho, L) :
    return rho

def normSchattenP(M, L, p) :
    dm = np.matrix(M, copy=True)
    for l in range(0, L) :
        dm[l, l] = 0.
    w, _ = np.linalg.eigh(dm)
    s = 0.
    for val in w :
        s += np.power(np.absolute(val), p)
    return np.power(s, 1./(p))

def normLP(M, L, p) :
    s = 0.
    for l1 in range(0, L) :
        for l2 in range(0, L) :
            if l1 != l2 :
                s += np.power(np.absolute(M[l1, l2]), p)
    return np.power(s, 1./(p))

def C2b(rho, x, L, N, p, norm, wrapper) :
    lhs = 0.
    for n in range(0, N) :
        pn = wrapper(PLL(x, rho, L, n), L)
        lhs += norm(pn, L, p)
    rhs = norm(wrapper(rho, L), L, p)
    return lhs <= rhs, lhs, rhs
