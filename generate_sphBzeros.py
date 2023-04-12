
import numpy as np
import fxstools.sphBzeros as sb
import sys
import os

#
# set the maximum l and number of zeros
#
lmax = 100
nt = 1000

#
# compute the zeros
#
z = sb.Jn_zeros(lmax, nt) 

outname = os.path.join( sys.path[0], "sphbzeros_lmax"+str(lmax)+"_nt"+str(nt)+".npy" )
np.save( outname, z )

print(z)
