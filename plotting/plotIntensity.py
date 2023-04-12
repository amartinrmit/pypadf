

import numpy as np
import matplotlib.pyplot as plt


#------------------------------------------------
#  INPUT PARAMETERS
#------------------------------------------------

#
# specify a path and filename
#
path = "/Users/andrewmartin/Work/Research/SF/codes/code_projects/pydiffraction/test/"
fname = "Intensity_water.txt"

fluence = 4
sigma = 15


#------------------------------------------------
#  PLOTTING CODE STARTS HERE
#------------------------------------------------
scale = 0.8490333   # this appears to be something to do with the non-linear time sampling; figured out "empirically"


# read in data
data = np.loadtxt( path+fname )


factor = fluence/(scale*sigma*np.sqrt(2*np.pi))
tsteps =  data[1:,1] - data[:-1,1]
print( "Total intensity :", np.sum(factor*data[0:-1,0]*tsteps))
print( "Total time :", np.sum(tsteps) )

print(data.shape)

plt.plot(tsteps)
plt.figure()

#
# plot the data
#
plt.plot( data[:,1], factor*data[:,0] )
plt.xlabel( "time (fs)" )
plt.ylabel( r'Beam Intensity (x10$^4$ J cm$^{-2}$ fs$^{-1}$)' )

plt.draw()
plt.show()
