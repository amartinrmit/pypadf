

import numpy as np
import matplotlib.pyplot as plt


#------------------------------------------------
#  INPUT PARAMETERS
#------------------------------------------------

#
# specify a path and filename
#
path = "/Users/andrewmartin/Work/Research/SF/codes/code_projects/pydiffraction/test/"
fname = "FormFactor_water2.txt"


#
# select some data to plot
#
# the times are in femtoseconds
# pelems contains the names of elements to plot
#
ptimes = [ -30, 10, 30]
pelems  = [ 'H', 'H', 'H'] 




#------------------------------------------------
#  PLOTTING CODE STARTS HERE
#------------------------------------------------

#
# read in the text file
#
input = np.loadtxt( open(path+fname, "r"), dtype="str" )
print( input.shape )
times = input[:,0].astype( np.float )
elements = input[:,1]
ffdata = input[:,2:].astype( np.float)

#
# count the unique elements
#
unique_elems = []
for e in elements:
    if e not in unique_elems:
        unique_elems.append( e )

npoints = int(input.shape[0]/len(unique_elems))   #number of time points

effdata = []
etimes = []
for i in range(len(unique_elems)):
    effdata.append( ffdata[i*npoints:(i+1)*npoints,:] )
    etimes.append( times[i*npoints:(i+1)*npoints] )
               #
# some information about the data
#
trange = np.max(times) - np.min(times)

print( "min and max times:", np.min(times), np.max(times) )
print( "number of time npoints:", int(input.shape[0]/2) )


#
# q data
#
nq = ffdata.shape[1]
q = np.arange( nq) * 0.1 * 2 * np.pi
invAngstromSymbol = r'$\mathrm{\AA}^{-1}$'



#
# plot all the data
#
plist = []
legend = []
for t, elem in zip( ptimes, pelems ):
    
    it = int( (t/trange)*npoints)
    if it<0:
        it = 0
    if it>npoints-1:
        it = npoints-1

    ie = unique_elems.index(elem)
    print( ie, elem, pelems )

    p, = plt.plot( q, effdata[ie][it,:] )

    legend += [ elem+" - t ="+str(t)+" fs" ]
    plist += [p]

plt.legend( plist, legend )
plt.xlabel( 'q ('+invAngstromSymbol+')' )
plt.ylabel( 'Atomic orm factor (average)' )

plt.draw()
plt.show()
 
