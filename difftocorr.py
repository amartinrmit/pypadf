
import numpy as np
import array
import os
import glob
import params.paramsCORR as params
import fxstools.padfio as padfio
import fxstools.correlation as correlation




if __name__ == '__main__': 
    

    print("\n-------------------------------------------------------------------------")
    print(" difftocorr.py : calculate correlation fucntion from diffraction patterns" )
    print("-------------------------------------------------------------------------")

    #
    # set up parameter class
    #
    p = params.paramsCORR()

    print("pypadf version ",p.version,"\n")

    #
    # Read input parameters from a file
    #
    p.read_config_file()
    p.qmax_calc()

    #
    # make the filelist
    #
    #(read in a file list or make a list using glob)
    p.load_flist_from_samplepath()

    #
    # read in a file and get the size of the file (nx,ny)
    #
    #print( "!DEBUG py3corrrelation nx ny", p.nx, p.ny)
    if p.nx==-1 or p.ny==-1:
        testimage = padfio.read_image(p.flist[0]) 
        p.nx, p.ny = testimage.shape[0], testimage.shape[1]
        #print( "!!!DEBUG py3corrrelation nx ny", p.nx, p.ny, testimage.shape)

    #
    # Set up an instance of the padf class
    #
    #print("DEBUG rebin:,", p.rebin)
    corr = correlation.correlation(path=p.outpath, tag=p.tag, flist=p.flist,
                 nx=p.nx, ny=p.ny, wl=p.wl, pw=p.pw, dz=p.dz, nth=p.nth,
                 nthreads=p.nthreads, npatterns=p.npatterns, 
                 bg_estimate=p.bg_estimate,
                 mask_flag=p.maskflag, crop_flag=p.cropflag, 
                 nxcrop=p.nxcrop, nycrop=p.nycrop,
                 dp_shift_flag=p.dp_shift_flag, 
                 shiftx=p.shiftx, shifty=p.shifty,
                 maskname=p.maskname, rebin=p.rebin, nstart=p.nstart,
                 diffcorr=p.diffcorrflag, outputdp=p.outputdp)
    #
    # calculate the padf
    #
    corrsum = np.zeros( (p.nx//2, p.nx//2, p.nth) )
    print("\nPerforming Correlations\n")
    corrsum = corr.calculate_correlation()

    print("\n")
    #
    # write the correlation function to file
    #
    #outname = p.outpath+p.tag+"_correlation_sum.dbin"  #append diff or bg as appropriate
    #padfio.write_dbin( outname, corrsum ) 
    outname = p.outpath / (p.tag+"_a_correlation_sum.npy")  #append diff or bg as appropriate
    np.save( outname, corrsum[:,:,:,0] ) 
    print("Written correlation sum:", outname)
    outname = p.outpath / (p.tag+"_b_correlation_sum.npy")  #append diff or bg as appropriate
    np.save( outname, corrsum[:,:,:,1] ) 
    print("Written correlation sum:", outname)
    print("Correlations Done.")
