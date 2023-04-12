

import numpy as np
import matplotlib.pyplot as plt


#------------------------------------------------
#  INPUT PARAMETERS
#------------------------------------------------

#
# specify a path and filename
#
path = "/Users/andrewmartin/Work/Research/SF/codes/code_projects/pydiffraction/test/"
fname = "MD_data.txt"




#------------------------------------------------
#  PLOTTING CODE STARTS HERE
#------------------------------------------------

# read in data
data = np.loadtxt( path+fname, skiprows=2 )

times = data[:,0]
electronE = data[:,-1]
electronN = data[:,-2]
charges = data[:,1:-2]

f = open( path+fname, "r" )
lines = f.readlines()
bits = lines[1].split()
zdata = np.array(bits[1:-2]).astype(np.int)
print(zdata)


#
# plot the trapped electron gas energy
#
plt.plot( times, electronE )
plt.xlabel( "time (fs)" )
plt.ylabel( "Electron Energy (eV)" )

#
# plot the trapped electron gas energy
#
invAsym3 = r'$\mathrm{\AA}^{3}$'
plt.figure()
plt.plot( times, electronN )
plt.xlabel( "time (fs)" )
plt.ylabel( r'Trapped Electron Density (e$^-$ per '+invAsym3+")" )


#
# plot the trapped electron gas energy
#
plt.figure()
plist = []
legend = []
for ie in range( charges.shape[1] ):
    p, = plt.plot( times, charges[:,ie] )
    plist.append( p )
    legend.append( "Z="+str(zdata[ie]) )

plt.legend( plist, legend )
plt.xlabel( "time (fs)" )
plt.ylabel( "Charge loss" )



plt.draw()
plt.show()
