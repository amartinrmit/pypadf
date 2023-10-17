

import numpy as np



diff = np.load('./demo/output/diff/hex_0.npy')

diff *= 1/diff

np.save('./demo/output/mask/hex_mask.npy', diff)

