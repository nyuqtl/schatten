import numpy as np
import qutip as qp

from simanneal import Annealer
from testing import normSchattenP, wrapperBiB, C2b, instanceX, testInstance

def rotate2D(z, a) :
    x = np.real(z)
    y = np.imag(z)
    x1 = x*np.cos(a) - y*np.sin(a)
    y1 = y*np.cos(a) + x*np.sin(a)
    return x1 + 1j*y1

def phase(inp, a) :
    inp = rotate2D(inp, a)
    m = np.linalg.norm(inp)
    inp = inp / m
    return inp

def amplitude(inp, coeff) :
    inp = inp*coeff
    m = np.linalg.norm(inp)
    inp = inp / m
    return inp

class OptimizeNorm(Annealer):
    def __init__(self, state, L, N, a, ma, testInstance=False):
        self.L = L
        self.N = N
        self.a = a
        self.ma = ma
        self.testInstance = testInstance
        super(OptimizeNorm, self).__init__(state)  # important!
    def move(self):
        """Randomly change one of measurements"""
        x, dm = self.state
        l = np.random.randint(0, self.L)
        n = np.random.randint(0, self.N)
        a = np.random.uniform(1e-8, self.a)
        ma = np.random.uniform(1e-8, self.ma)
        x[l,:] = phase(x[l,:], a)
        x[l,:] = amplitude(x[l,:], ma)
        if self.testInstance :
            if not testInstance(x, self.L, self.N) :
                raise ValueError('Invalid measurement kets.')
        #dm = np.matrix(qp.rand_dm(self.L).full())
    def energy(self):
        """Calculates the difference between RHS and LHS of C2(b)"""
        x, dm = self.state
        _, lhs, rhs = C2b(dm, x, self.L, self.N, 1., normSchattenP, wrapperBiB)
        return rhs - lhs

# choose dimensions
L = 4
N = 7

# L is dimension of density matrix
# N is number of measurements

# measurements can be represented by L kets of size N
# (N outcomes)

# get some random density matrix
dm = np.matrix(qp.rand_dm(L).full())

# get set of random kets that will be used for measurement
x = instanceX(L, N)

print 'random sampling for density matrix for which C2b vanishes\n'
_, lhs, rhs = C2b(dm, x, L, N, 1., normSchattenP, wrapperBiB)
mini = rhs - lhs
for i in range(0, 10000) :
    rho = np.matrix(qp.rand_dm(L).full())
    _, lhs, rhs = C2b(rho, x, L, N, 1., normSchattenP, wrapperBiB)
    nmini = rhs - lhs
    if nmini < mini :
        print nmini
        mini = nmini
        dm = rho
        if np.isclose([mini], [0.])[0] :
            print '\nfound tightest C2(b) bound which is approximately 0.'
            break

print 'beginning simulated annealing to find measurement violating C2(b) as local optimum\n'
initial_state = tuple([x, dm])
opt = OptimizeNorm(initial_state, L, N, .5*np.pi, .5, testInstance=False)
opt.Tmax = 200000.0
final_state, difference = opt.anneal()

print '\n\n'

if np.isclose([difference], [0.])[0] :
    print 'no C2(b) violation found'
