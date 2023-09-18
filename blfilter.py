
import numpy as np
import array
import os
import params.paramsFILT as params
import fxstools.padflib as padflib
import fxstools.padfio as padfio


print("\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!")
print(" blfilter.py : blur a model volume to match experimental resolution")
print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!")


#
# set up parameter class
#
p = params.paramsFILT()


#
# Read input parameters from a file
#
p.read_config_file()
p.readCorr = False

#
# Set up an instance of the padf class
#
#print("DEBUG pyblfilter.py qmax, rmax", p.qmax, p.rmax)
#padf = pld.padfcl( nl=p.nl, nr=p.nr, nq=p.nq, qmin=p.qmin, qmax=p.qmax, rmax=p.rmax  )

fname = p.path_to_string( p.padffile )
padf = padfio.read_image( fname )


pcl = padflib.padfcl(p.nl+2,p.nlmin,p.nr,p.nth,p.nq,p.qmin,p.qmax,p.rmax,padf,wl=p.wl) 


nl,nlmin,nr,nth,nq,qmin,qmax,rmax,wl = p.nl+2,p.nlmin,p.nr,p.nth,p.nq,p.qmin,p.qmax,p.rmax,p.wl

print("Min and max q values:",p.qmin, p.qmax)
print("Min and max l values", p.nlmin, nl)

#
# calculate Blrr
#
blrr = pcl.calcBlrr( padf, nl )

#
# use dsB matrix to map to Fourier space
#
print("Calculating blqq matrices")
dsBlist = []
dsBinvlist = []
#print("nlmin, nl", nlmin, nl) #DEBUG
for l in range(nlmin, nl, 1):
    nmax = pcl.blqq.sphB_samp_nmax( l, rmax, qmax)
    #print("TEST!!!!! l nmax", l, nmax)  #DEBUG
    dsBlist.append( pcl.dsBmatrix( l, nmax ) )
    dsBinvlist.append( pcl.dsBmatrix_inv( l, nmax, qmin ) )

for l in range( nlmin, nl, 1 ):
    pcl.blqq.l[l].data = np.dot(np.dot( dsBinvlist[l], blrr[l] ), dsBinvlist[l].transpose() )

pcl.blqq.l[0].data *= 0.0 #to ensure l=0 is removed


## apply the Ewald filter
####blqqfilt = pcl.Blqq_ewald_filter( pcl.blqq, qwid=qwid )
####blqqfilt = pcl.blqq

if p.blqqtocorr:
    ### calculate a correlatino function from the Blqq matrices
    print("Converting blqq to correlation function, then recalculate blqq")
    if p.readCorr:
        pcl.corrvol = corrin
    else:
        corr = pcl.Blqq_to_corr_fast( pcl.blqq, interpolate=p.interpolate, order=p.order)
        """
        disp = np.zeros( (corr.shape[0], corr.shape[2] ) )
        for i in np.arange(corr.shape[0]):
            disp[i,:] = corr[i,i,:]
        plt.imshow( disp, origin='lower')
        plt.show()
        exit()
        """
        pcl.corrvol = corr #corrin
    if p.interpolate:
        pcl.Blqq_calc_fast_interp(order=p.order)
    else:
        pcl.Blqq_calc_fast(outputfmat=True)
blqqfilt = pcl.blqq

#
# invert dsB matrix
#
print("Converting Blqq to Blrr")
for l in range( nlmin, nl, 1 ):
    blrr[l] = np.dot(np.dot( dsBlist[l], blqqfilt.l[l].data ), dsBlist[l].transpose() )

#
# map Blrr to PADF
#
print("Converting Blrr into PADF")
padf_filt = pcl.Blrr_to_padf( blrr, padf.shape )


print("Calculation finished!")
#
# output the filtered PADF
#
print( "PADF filtered min & max:", padf_filt.min(), padf_filt.max() )
outname = p.path_to_string(p.outpath / (p.tag+"_padf_filt.npy"))
np.save( outname, padf_filt )

outname = p.path_to_string(p.outpath / (p.tag+"_corr_predict.npy"))
np.save( outname, pcl.corrvol )

print("\n Output file:")
print(outname)








"""
# *** OLD CODE ***

#
#  Calculate the filter files
#
print(p.nl, "- nl")
print(" ")
print(" ")
dsBlist = []
dsBinvlist = []
for l in range(p.nlmin, p.nl, 1):

    nmax = padf.blqq.sphB_samp_nmax( l, p.rmax, p.qmax)
    print("l nmax", l, nmax)

    dsB = padf.dsBmatrix( l, nmax )
    dsBinv = padf.dsBmatrix_inv( l, nmax, p.qmin )

    filtermat = np.dot( dsB, dsBinv)

    dsBlist.append(dsB)
    dsBinvlist.append(dsBinv)

    #
    # output filter files  
    #
    if p.savefilt:
        outname = p.outpath+p.tag+"_l"+str(l)+"_filter.npy"
        np.save( outname, filtermat )

        outname = p.outpath+p.tag+"_l"+str(l)+"_dsB.npy"
        np.save( outname, dsB )

#
# filter the padf 
#

#
# read the correlation file
#
if p.filterpadf:
    if p.padffile[-4:]=='dbin':
        padfin = padfio.read_correlation( p.padffile, 0 )
        padfin = padfin.reshape( p.nq, p.nq, p.nth )
    elif p.padffile[-3:]=='npy':
        padfin = np.load(p.padffile)
    else:
        print("Correlation volume must be dbin or npy format", p.padffile)
        exit()


    # calculate Blrr
    blrr = padf.calcBlrr( padfin, p.nl )

    
    for l in range( p.nlmin, p.nl, 1 ):
        padf.blqq.l[l].data = np.dot(np.dot( dsBinvlist[l], blrr[l] ), dsBinvlist[l].transpose() )

    if p.blqqtocorr:
        padf.corrvol = padf.Blqq_to_corr_fast( padf.blqq )
        padf.Blqq_calc_fast(outputfmat=True)


    # invert dsB matrix
    for l in range( p.nlmin, p.nl, 1 ):
        blrr[l] = np.dot(np.dot( dsBlist[l], padf.blqq.l[l].data ), dsBlist[l].transpose() )

    padf.padf = padf.Blrr_to_padf( blrr, padf.shape )


    # output fliltered padf
    outname = p.outpath+p.tag+"_padf.dbin"
    padfio.write_dbin( outname, padf.padf )
    outname = p.outpath+p.tag+"_padf.npy"
    np.save(outname, padf.padf)
"""
