
import numpy as np
#import array
import os
import params.paramsPADF as params
import fxstools.padflib as plib
import fxstools.padfio as padfio




#
# set up parameter class
#
p = params.paramsPADF()

print("\n")
print(" ____   __   ____  ____")
print("(  _ \\ / _\\ (    \\(  __)")
print(" ) __//    \\ ) D ( ) _) ")
print("(__)  \\_/\\_/(____/(__) ")
print("\n")
print("pypadf version ",p.version,"\n")

#
# Read input parameters from a file
#
p.read_config_file()

#
# read the correlation file
#
#help(padfio)
corrfile = p.path_to_string(p.corrfile)
if corrfile[-4:]=='dbin':
    corrvol = padfio.read_correlation( corrfile, 0 )
    corrvol = corrvol.reshape( p.nq, p.nq, p.nth )
elif corrfile[-3:]=='npy':
    corrvol = np.load(corrfile)
else:
    print("Correlation volume must be dbin or npy format", corrfile)
    exit()

#
# normalise by beamnorm
#
corrvol *= 1.0/(p.beamnorm**2)



#
# Set up an instance of the padf class
#
# (2 is added to nl so that nl represents 2x the actual number of spherical harmonics used)
#
padf = plib.padfcl( nl=p.nl+2, nlmin=p.nlmin, nr=p.nr, nq=p.nq, 
                   qmin=p.qmin, qmax=p.qmax, rmax=p.rmax, nth=p.nth,
                   corrvol=corrvol, wl=p.wl, method=p.method,legendre_norm=p.legendre_norm)

#
# calculate the padf
#
padf.padf_calc()

#
# change units from m^-6 to Angstroms^-6
#
padf.padf *= 10e-60


#
# normalise by density
#
padf.normalise_padf_with_density( p.density )


#
# write the padf to file
#
#outname = p.outpath+p.tag+"_padf.dbin"
#padfio.write_dbin( outname, padf.padf )
outname = p.path_to_string(p.outpath / (p.tag+"_padf.npy"))
np.save(outname, padf.padf)
#
# write out the blqq/blrr arrays for debugging
#
blpath = p.outpath / "blmatrices"
blstr = p.path_to_string(blpath)
if os.path.exists(blstr)!=True:
    os.mkdir(blstr)
padf.blqq.write_blqq_array_as_npy( blpath, p.tag+"_blqq" )
#padf.blqqtmp.write_blqq_array_as_npy( p.outpath, p.tag+"_blqqtmp" )
padf.blrr.write_blqq_array_as_npy( blpath, p.tag+"_blrr" )
