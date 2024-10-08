

import numpy as np
import matplotlib.pyplot as plt
import os
import time

import params.paramsDIFF_batch as params
import fxstools.diffraction as df
import fxstools.quaternions as quaternions
import fxstools.pydiffractionio as pio

#
# set up parameter class and read all parameters from file
#
p = params.paramsDIFFBATCH()

print("saxs.py")
print("pypadf version ",p.version,"\n")

#
# Read parameters from the configuration file
#
p.read_diff_parameters_from_file()


#
# Set up the diffraction class object
#
pdbname = str(p.pdbname.resolve())
#outpath = str(p.outpath.resolve())
difrct = df.diffraction(pdbname, p.outpath, p.tag, p.fext, p.nx, p.wl, p.dz, p.pw, p.cenx,
                        p.ceny, p.henkeflag, np.array([p.rx,p.ry,p.rz]),p.rtheta,p.rotflag,
                        npulse=p.npulse,beamarea=p.beamarea )

with open(p.outname,'a') as f:
    f.write(f"qmax = {difrct.q1d[-1]}\n")
    f.write(f"elements =")
    for ie in difrct.pdb.elements:
        f.write(f' {ie}')

difrct.pdb.maxdims()
mins =  (difrct.pdb.xmin,difrct.pdb.ymin,difrct.pdb.zmin)

if(p.alen<difrct.pdb.xlen)or(p.blen<difrct.pdb.ylen)or(p.clen<difrct.pdb.zlen):
    uc = df.unitcell(a=difrct.pdb.xlen,b=difrct.pdb.ylen,c=difrct.pdb.zlen)
else:
    uc = df.unitcell(a=p.alen,b=p.blen,c=p.clen)

#save the pdbname
pdbsave = pdbname
 
#
# calculate the saxs pattern
#
difrct.load_pdb()

start = time.perf_counter()
difrct.saxs(nr=p.saxs_nr, box_nx=p.saxs_box_nx, norm_rmin=p.saxs_norm_rmin, norm_rmax=p.saxs_norm_rmax)
end = time.perf_counter()
print("Time to calculate the saxs pattern:", end-start, "seconds")
fname = difrct.outpath / (difrct.tag+"_saxs"+difrct.fext)
pio.write_image( str(fname.resolve()), difrct.saxs) 

fname = difrct.outpath / (difrct.tag+"_saxs_partial"+difrct.fext)
pio.write_image( str(fname.resolve()), difrct.saxs_partial) 

fname = difrct.outpath / (difrct.tag+"_pairdist"+difrct.fext)
pio.write_image( str(fname.resolve()), difrct.partial_pd) 

fname = difrct.outpath / (difrct.tag+"_boxpd"+difrct.fext)
pio.write_image( str(fname.resolve()), difrct.pdbox) 


#
# display the saxs pattern
#
if p.display==True:
    #
    # Plot the output to screen
    #

    # q data
    nq = p.nx
    q = np.arange( nq) * 0.1 * 2 * np.pi
    invAngstromSymbol = r'$\mathrm{\AA}^{-1}$'


    gamma = 1.0
    plt.plot( difrct.q1d, difrct.saxs**gamma )

    plt.xlabel( 'q ('+invAngstromSymbol+')' )
    plt.ylabel( 'Intensity (arb. units)' )
    plt.ylim([0,np.max(difrct.dp1d)**gamma])
    plt.xlim([0,np.max(difrct.q1d)])

    plt.draw()
    plt.show()

