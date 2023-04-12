"""

Plotting fluctuation scattering volumes

q-space (intensity correlation function) or real-space (PADF)

"""

import numpy as np
import os
import matplotlib.pyplot as plt
import params.paramsPLOT as params
import fxstools.padfplot as pp


#
# set up parameter class
#
p = params.paramsPLOT()

print("pypadf version ",p.version,"\n")

#
# Read input parameters from a file
#
p.read_config_file()

#
# set up plotting dimensions
#
pdims = pp.padfplot_dims(rmin=p.rmin,rmax=p.rmax,rval=p.rval,rval2=p.rpval,thval=p.thval,rwid=p.rwid)

#
# read the correlation file
#
#help(padfio)
fname = p.path_to_string(p.fname)
if fname[-4:]=='dbin':
    fxsvol = padfio.read_correlation( fname, 0 )
    fxsvol = corrvol.reshape( p.nq, p.nq, p.nth )
elif fname[-3:]=='npy':
    fxsvol = np.load(fname)
else:
    print("Correlation volume must be dbin or npy format", p.fname)


# scale radial dimension
if p.power!=0:
    fxsvol = pp.mult_radial_polynomial( fxsvol, p.power, p.rmin, p.rmax )

# subtract angular mean 
if p.submean==True:
    fxsvol = pp.remove_angular_average( fxsvol )

# multiply by  |sin(theta)|
if p.sintheta==True:
    fxsvol = pp.multiply_by_sintheta( fxsvol, p.thmin, p.thmax )

# convolve the volume
if p.convolve==True:
    cr = pdims.get_ir(fxsvol.shape[0], p.rwid, flt=True)
    cth = pdims.get_ith(fxsvol.shape[2],p.thwid,flt=True)
    fxsvol = pp.convolve_gaussian_n( fxsvol, rad=cr, rady=cr, radz=cth)

#extract the slice or line
disp = pp.extract_section( fxsvol, pdims, p.stype ) 

#
# Plot the section of the 3D volume
#
rlabel, thlabel = pp.generate_unit_labels(p.runits,p.rq)
r1 = pdims.get_ir(fxsvol.shape[0], p.rmaxdisp)
r0 = pdims.get_ir(fxsvol.shape[0], p.rmindisp)

if p.stype=='reqr' or p.stype=='rconst':
    plt.imshow(disp[r0:r1,:], origin='lower', extent=[0,p.thmax,p.rmindisp,p.rmaxdisp], aspect=360/(p.rmaxdisp-p.rmindisp), interpolation='gaussian')
    if (p.chigh>0) and (p.clow>=0):
        plt.clim([p.clow, p.chigh])
    else:
        plt.clim([np.min(disp[r0:r1,:])*p.scalel,np.max(disp[r0:r1,:])*p.scale])
    plt.colorbar()
    plt.ylabel(rlabel)
    plt.xlabel(thlabel)

elif p.stype=='thconst':
    plt.imshow(disp[r0:r1,r0:r1], origin='lower', extent=[p.rmindisp,p.rmaxdisp,p.rmindisp,p.rmaxdisp], aspect=1, interpolation='gaussian')
    if (p.chigh>0) and (p.clow>=0):
        plt.clim([p.clow, p.chigh])
    else:
        plt.clim([np.min(disp[r0:r1,r0:r1])*p.scalel,np.max(disp[r0:r1,r0:r1])*p.scale])
    plt.colorbar()
    plt.ylabel(rlabel)
    plt.xlabel(rlabel)
    
elif p.stype=="rline":
    rvals = (p.rmaxdisp-p.rmindisp)*np.arange(r1-r0)/(r1-r0) + p.rmindisp
    plt.plot( rvals, disp[r0:r1] )
    if p.rq=='r':
        plt.ylabel("PADF (arbitrary units)")
    elif p.rq=='q':
        plt.ylabel("Intensity correlation function (arb. units)") 
    plt.xlabel(rlabel)

elif p.stype=="thline":
    thvals = (p.thmax/disp.size)*np.arange(disp.size)
    plt.plot( thvals, disp )
    if p.rq=='r':
        plt.ylabel("PADF (arbitrary units)")
    elif p.rq=='q':
        plt.ylabel("Intensity correlation function (arb. units)") 
    plt.xlabel(thlabel)

# save the figure data
outname = p.path_to_string( p.outpath / (p.fname_basenoext+"_"+p.suffix+"_"+p.stype+".npy") )
np.save(outname, disp)
# save the figure image
outname = p.path_to_string( p.outpath / (p.fname_basenoext+"_"+p.suffix+"_"+p.stype+".png") )
plt.savefig(outname)

plt.draw()
plt.show()
