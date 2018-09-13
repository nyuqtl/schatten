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

def testKets(kets) :
    for ket in kets :
        s = 0.
        for k in ket :
            s += np.power(np.abs(k), 2.)
        if not np.isclose([s], [1.])[0] :
            return False
    return True

def dm(kets, prob, M) :
    dm = prob[0] * qp.ket2dm(kets[0])
    for i in range(1, M) :
        dm += prob[i] * qp.ket2dm(kets[i])
    return np.matrix(dm.full())

def randomChange(kets, prob, N, M, a, ma, mx) :
    # pick random ket
    index1 = np.random.randint(M)
    # pick random element in this ket
    index2 = np.random.randint(N)
    # change phase
    kets[index1] = phase(kets[index1], index2, a)
    if not testKets(kets) :
        raise Exception('phase made wrong ket')
    # change amplitude
    kets[index1] = amplitude(kets[index1], index2, ma)
    if not testKets(kets) :
        raise Exception('amplitude made wrong ket')
    # change probabilities
    ch = np.arange(M)
    np.random.shuffle(ch)
    prob[ch[0]], prob[ch[1]] = moveProb(prob[ch[0]], prob[ch[1]], mx)
    return kets, prob

def rotate2D(z, a) :
    x = np.real(z)
    y = np.imag(z)
    x1 = x*np.cos(a) - y*np.sin(a)
    y1 = y*np.cos(a) + x*np.sin(a)
    return x1 + 1j*y1

def phase(ket, index, a) :
    inp = ket.full()
    inp[index] = rotate2D(inp[index], a)
    m = np.linalg.norm(inp)
    inp = inp / m
    return qp.Qobj(inp)

def amplitude(ket, index, coeff) :
    inp = ket.full()
    inp[index] = ket[index]*coeff
    m = np.linalg.norm(inp)
    inp = inp / m
    return qp.Qobj(inp)

def moveProb(p1, p2, mx) :
    r = np.random.uniform(low=0.0, high=mx, size=(1,))[0]
    return p1 - r, p2 + r

N = 4
M = 5

# make random ket
kets = []
prob = np.random.dirichlet(np.ones(M))
for i in range(0, 5) :
    ket = qp.rand_ket(N)
    kets.append(ket)


print np.sum(prob)
print testKets(kets)
rho = dm(kets, prob, M)
testDM(rho, silent=False)

kets2, prob2 = randomChange(kets, prob, N, M, 1.01, np.pi, 0.2)
rho2 = dm(kets2, prob2, M)
print '---'
print np.sum(prob2)
print testKets(kets2)
testDM(rho2, silent=False)





'''
ket = qp.rand_ket(2).full()
z1 = ket[0][0]
z2 = ket[1][0]
print z1, z2
print np.power(np.absolute(z1), 2.) + np.power(np.absolute(z2), 2.)
nz1, nz2 = moveComplex(z1, z2, 0.12, 0.13)
print nz1, nz2
print np.power(np.absolute(nz1), 2.) + np.power(np.absolute(nz2), 2.)
'''

#randomChange(kets, prob)

'''
ket = qp.rand_ket(5)
print ket.full()
ket = phase(ket, 0, 1.2)
ket = amplitude(ket, 0, 1.1)
print ket.norm()
print qp.Qobj(ket).full()
'''

# small random change in the ket
#H = qp.rand_herm(N)

# make random hermitian operator


# small random change in the hermitian operator
