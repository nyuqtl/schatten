import cPickle as pickle
import io

import numpy as np
import qutip as qp

from tqdm import tqdm
from tools import instanceX

m = 1000

for L in [2, 4, 8] :
    for N in [2, 4, 8] :
        print 'L = %s, N = %s' % (str(L), str(N))
        filename = 'instances/x_%s_%s_%s.p' % (str(L), str(N), str(m))
        with io.open(filename, 'wb') as f :
            pickler = pickle.Pickler(f)
            for l in tqdm(range(m)) :
                x = instanceX(L, N)
                pickler.dump(x)
