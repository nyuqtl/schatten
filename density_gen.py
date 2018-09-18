import cPickle as pickle
import io

import numpy as np
import qutip as qp

from tqdm import tqdm

m = 1000

for L in [2, 4, 8] :
    print 'L = %s' % str(L)
    filename = 'instances/rho_%s_%s.p' % (str(L), str(m))
    with io.open(filename, 'wb') as f :
        pickler = pickle.Pickler(f)
        for l in tqdm(range(m)) :
            rho = np.matrix(qp.rand_dm(L).full())
            pickler.dump(rho)
