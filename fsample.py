import cPickle as pickle
import io
import numpy as np

from tqdm import tqdm
from testing import C2b, wrapperBiB, wrapperIdentity, normLP, normSchattenP

m = 1000

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
for L in [2] :
    for N in [2] :
        print 'L = %s, N = %s' % (str(L), str(N))

        file_rhos = 'instances/rho_%s_%s.p' % (str(L), str(m))
        file_xs = 'instances/x_%s_%s_%s.p' % (str(L), str(N), str(m))

        rhos = io.open(file_rhos, 'rb')
        xs = io.open(file_xs, 'rb')

        unpickler_rhos = pickle.Unpickler(rhos)
        unpickler_xs  = pickle.Unpickler(xs)

        for l in tqdm(range(m)) :
            rho = unpickler_rhos.load()
            x = unpickler_xs.load()
            LP_ide, _, _ = C2b(rho, x, L, N, p, normLP, wrapperIdentity)
            LP_bib, _, _ = C2b(rho, x, L, N, p, normLP, wrapperBiB)
            sc_ide, _, _ = C2b(rho, x, L, N, p, normSchattenP, wrapperIdentity)
            sc_bib, _, _ = C2b(rho, x, L, N, p, normSchattenP, wrapperBiB)
            counter['LP']['ide'] += int(not LP_ide)
            counter['LP']['bib'] += int(not LP_bib)
            counter['sc']['ide'] += int(not sc_ide)
            counter['sc']['bib'] += int(not sc_bib)

        rhos.close()
        xs.close()

line1 = ('%02d' % int(p), '{:11d}'.format(counter['LP']['ide']), '{:12d}'.format(counter['LP']['bib']))
line2 = ('%02d' % int(p), '{:11d}'.format(counter['sc']['ide']), '{:12d}'.format(counter['sc']['bib']))

print '\nnumber of violating instances:\n'
print '+--- norm ----+-- Pn / P --+-- B + iBd --+'
print '| L%s         |%s |%s |' % line1
print '| Schatten-%s |%s |%s |' % line2
print '+-------------+------------+-------------+\n'
