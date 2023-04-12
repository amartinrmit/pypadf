

import numpy as np
import matplotlib.pyplot as plt


#
# paths containing the data
#
paths = [ "replace with path1",
         "replace with path2",
         "etc." ]

tags =  [ "replace with tag1",
         "replace with tag2",
         "etc." ]

gamma = 0.3

refname = path[0]+tag[0]+"1D_diffraction.txt"
refdata = np.loadtxt( refname )

for path, tag in zip( paths, tags ):

    fname = path+tag+"1D_diffraction.txt"
    data = np.loadtxt( fname )

	
    # plot the data as normal
    plt.plot( data[0,:], data[1,:]**0.3 )

    #

invAngstromSymbol = r'$\mathrm{\AA}^{-1}$'
plt.xlabel( 'q ('+invAngstromSymbol+')' )
plt.ylabel( 'Intensity (arb. units)' )
plt.ylim([0,np.max(data[1,:])**gamma])
plt.xlim([0,np.max(data[0,:])])

plt.draw()
plt.show()

    
