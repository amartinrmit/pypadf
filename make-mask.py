

import numpy as np



diff = np.load('./tmp/diff/hex_0.npy')

diff *= 1/diff

np.save('./tmp/mask/hex_mask.npy', diff)

