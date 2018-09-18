import qutip as qp
import numpy as np

from tqdm import tqdm
from testing import C2b, wrapperBiB, wrapperIdentity, normLP, normSchattenP
from tools import instanceX

total = 0
all = 0

counter = {
    'LP': {
        'ide': 0,
        'bib': 0
    },
    'sc': {
        'ide': 0,
        'bib': 0
    }
}

p = 1.
m = 100000

lvals = [2, 4, 8, 16, 32]
nvals = [2, 4]

for L in lvals :
    for N in nvals :
        print 'L = %s, N = %s' % (str(L), str(N))
        for l in tqdm(range(m)) :
            rho = np.matrix(qp.rand_dm(L).full())
            x = instanceX(L, N)
            LP_ide, _, _ = C2b(rho, x, L, N, p, normLP, wrapperIdentity)
            LP_bib, _, _ = C2b(rho, x, L, N, p, normLP, wrapperBiB)
            sc_ide, _, _ = C2b(rho, x, L, N, p, normSchattenP, wrapperIdentity)
            sc_bib, _, _ = C2b(rho, x, L, N, p, normSchattenP, wrapperBiB)
            counter['LP']['ide'] += int(not LP_ide)
            counter['LP']['bib'] += int(not LP_bib)
            counter['sc']['ide'] += int(not sc_ide)
            counter['sc']['bib'] += int(not sc_bib)
            total += int(not LP_ide)
            total += int(not LP_bib)
            total += int(not sc_ide)
            total += int(not sc_bib)
        all += m
        print '\ntotal violations so far: %s\n' % str(total)

line1 = ('%02d' % int(p), '{:11d}'.format(counter['LP']['ide']), '{:12d}'.format(counter['LP']['bib']))
line2 = ('%02d' % int(p), '{:11d}'.format(counter['sc']['ide']), '{:12d}'.format(counter['sc']['bib']))

print '\ntotal number of instances: %s\n' % str(all)
print '\ntotal number of violated instances: %s\n' % str(total)
print '\nnumber of violating instances per test:\n'
print '+--- norm ----+-- Pn / P --+-- B + iBd --+'
print '| L%s         |%s |%s |' % line1
print '| Schatten-%s |%s |%s |' % line2
print '+-------------+------------+-------------+\n'
