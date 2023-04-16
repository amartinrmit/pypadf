

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

print("diffract.py")
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

difrct.pdb.maxdims()
mins =  (difrct.pdb.xmin,difrct.pdb.ymin,difrct.pdb.zmin)

if(p.alen<difrct.pdb.xlen)or(p.blen<difrct.pdb.ylen)or(p.blen<difrct.pdb.zlen):
    uc = df.unitcell(a=difrct.pdb.xlen,b=difrct.pdb.ylen,c=difrct.pdb.zlen)
else:
    uc = df.unitcell(a=p.alen,b=p.blen,c=p.clen)

#save the pdbname
pdbsave = pdbname
 
pdbname  = pdbname[:-4]+"_shifted.pdb"

#
# calculate the 2D diffraction pattern
#
for i in np.arange(p.npatterns):
    difrct.axis, difrct.theta = quaternions.random_vec_angle() 

    # shift the pdb coordinates using a unit cell
    difrct.circ_shift_coordinates( uc, mins,True)
    # write a new pdb file with the shifted coordinates
    #outname = p.pdbname[:-4]+"_shifted.pdb"
    difrct.pdb.write_pdb(pdbname)
    # load the pdb file, to ensure new atom infomration is appropriately sorted
    difrct.load_pdb()

    start = time.perf_counter()
    difrct.diffraction2D()
    if p.polarisation is True: difrct.dp2d *= difrct.polarisation_factor([p.px, p.py])
    if p.poisson is True: difrct.dp2d = difrct.poisson_sample(difrct.dp2d)
    end = time.perf_counter()
    print( i+1,"/",p.npatterns,"  Time to calculate 2D diffraction pattern:", end-start, "seconds", end="\r")

    #
    # output 2D diffraction pattern to file
    #
    fname = difrct.outpath / (difrct.tag+"_"+str(i)+difrct.fext)
    #print("Output file:", fname)
    pio.write_image( str(fname.resolve()), difrct.dp2d )

    #
    # Calculate 1D pattern
    #
    difrct.diffraction1D()



#
# Plot the output to screen
#

# q data
nq = p.nx
q = np.arange( nq) * 0.1 * 2 * np.pi
invAngstromSymbol = r'$\mathrm{\AA}^{-1}$'


gamma = 1.0
#plt.imshow( np.log(difrct.dp2d) )
plt.imshow( np.abs(difrct.dp2d)**gamma )
#plt.imshow( difrct.sflist[2].sf2d)


plt.figure()
plt.plot( difrct.q1d, difrct.dp1d**gamma )

plt.xlabel( 'q ('+invAngstromSymbol+')' )
plt.ylabel( 'Intensity (arb. units)' )
plt.ylim([0,np.max(difrct.dp1d)**gamma])
plt.xlim([0,np.max(difrct.q1d)])

plt.draw()
plt.show()

