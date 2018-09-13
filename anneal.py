import numpy as np
import qutip as qp

def testDM(rho, silent=True) :
    w, _ = np.linalg.eigh(rho)
    trace = np.trace(rho)
    t1 = np.isclose([trace], [1.])[0]
    t2 = (w >= 0.).all()
    t3 = (rho.H == rho).all()
    if not silent :
        print 'trace is one: %s' % str(t1)
        print 'eigenvalues positive: %s' % str(t2)
        print 'is Hermitian: %s' % str(t3)
    return t1 and t2 and t3

def testKets(kets) :
    for ket in kets :
        s = 0.
        for k in ket :
            s += np.power(np.abs(k), 2.)
        if not np.isclose([s], [1.])[0] :
            return False
    return True

def dm(kets, prob, N) :
    dm = qp.Qobj(np.zeros((N, N), np.float64))
    for ket, p in zip(kets, prob) :
        dm += p * (ket * ket.dag())
    return np.matrix(dm.full())

#def randomChange(kets, prob) :

def rotate2D(x, y, a) :
    x1 = x*np.cos(a) - y*np.sin(a)
    y1 = y*np.cos(a) + x*np.sin(a)
    return x1, y1

def moveComplex(z1, z2, ma, mx) :
    a1 = np.random.uniform(low=0.0, high=ma, size=(1,))[0]
    a2 = np.random.uniform(low=0.0, high=ma, size=(1,))[0]
    r = np.random.uniform(low=0.0, high=mx, size=(1,))[0]
    x1, y1 = rotate2D(np.real(z1), np.imag(z1), a1)
    x2, y2 = rotate2D(np.real(z2), np.imag(z2), a2)
    m1 = np.sqrt(np.power(x1, 2.) + np.power(y1, 2.))
    m2 = np.sqrt(np.power(x2, 2.) + np.power(y2, 2.))
    nz1 = (x1/m1)*(m1+r) + 1j*(y1/m1)*(m1+r)
    nz2 = (x2/m2)*(m2-r) + 1j*(y2/m2)*(m2+r)
    return nz1, nz2

def moveProb(p1, p2, mx) :
    r = np.random.uniform(low=0.0, high=mx, size=(1,))[0]
    return p1 - r, p2 + r

N = 4
M = 5

# make random ket
kets = []
prob = np.random.dirichlet(np.ones(N))

for i in range(0, 5) :
    ket = qp.rand_ket(N)
    kets.append(ket)

rho = dm(kets, prob, N)

#testDM(rho, silent=False)

print testKets(kets)


# small random change in the ket
#H = qp.rand_herm(N)

# make random hermitian operator


# small random change in the hermitian operator
