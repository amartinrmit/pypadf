"""

display a 2D diffraction pattern

"""

import numpy as np
import matplotlib.pyplot as plt
import params.params as params
import fxstools.padfio as io
import os

# scale a function by gamma
def corr_rescale( plane, gamma ):
     disp = plane*0.0
     ihigh = np.where( plane > 0 )
     ilow = np.where( plane < 0 )
     disp[ihigh] = np.abs( plane[ihigh] )**gamma
     disp[ilow] = - np.abs( plane[ilow] )**gamma
     return disp

#
# Make an instance of the parameter class an add key parameters
#
p = params.params()
ch = ["PLOTDIFFRACTION"]
p.add_parameter("fname", "str", cmdline="-f", cmdline2="--fname",help="Name of the diffraction file",nargs=1,header="DEFAULT",pathflag=False)
p.add_parameter("scale", 1.0, cmdline="--scale",cmdline2="-sc",help="upper clim val in np.max(image)*scale", nargs=1,header=ch[0],pathflag=False)
p.add_parameter("scalel", 1.0, cmdline="--scalel",cmdline2="-scl", help="lower clim val in np.min(image)*scalel", nargs=1,header=ch[0],pathflag=False)
p.add_parameter("clow", -1.0, cmdline="--clow",cmdline2="-cl", help="absolute clim lower limit (priotity over scale)", nargs=1,header=ch[0],pathflag=False)
p.add_parameter("chigh", -1.0, cmdline="--chigh",cmdline2="-ch", help="absolute clim upper limit (priority over scalel)", nargs=1,header=ch[0],pathflag=False)
p.add_parameter("log", False, cmdline="--log",cmdline2="-l", help="log scale on all plots", nargs=1, header=ch[0],pathflag=False)
p.add_parameter("gamma", -1.0, cmdline="--gamma",cmdline2="-g", help="gamma value to apply to all images", nargs=1, header=ch[0],pathflag=False)


#
# read the config file (and command line parameters)
#
#args = p.parser.parse_args()
#print("Test command line reading:", args.fname[0], args.scale)
p.parse_commandline_args()
#help(p.fname)


#
# read in the diffraction pattern
#
data = io.read_image(p.fname)

#
# scale the image by gamma value if required
#
if p.gamma > 0.0:
    data = corr_rescale(data, p.gamma)    

#
# plot the 2D function
#
plt.imshow(data, origin='lower')
if (p.chigh>0) and (p.clow>=0):
    plt.clim([p.clow, p.chigh])
else:
    plt.clim([np.min(data)*p.scalel,np.max(data)*p.scale])
plt.colorbar()
#plt.ylabel(rlabel)

plt.draw()
plt.show()


