

import numpy as np
import sys



diff = np.load(sys.argv[1])

diff *= 1/diff

np.save(sys.argv[2], diff)

