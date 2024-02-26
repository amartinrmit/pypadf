
import numpy as np
import array
import os
import glob
import params.paramsCORR as params
import fxstools.padfio as padfio
import fxstools.correlation as correlation




if __name__ == '__main__':

    print("\n-------------------------------------------------------------------------")
    print(" difftocorr.py : calculate correlation function from diffraction patterns" )
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
    # Set up an instance of the correlation class
    #
    #print("DEBUG rebin:,", p.rebin)
    corr = correlation.correlation(path=p.outpath, tag=p.tag, flist=p.flist,
                 nx=p.nx, ny=p.ny, wl=p.wl, pw=p.pw, dz=p.dz, nth=p.nth,
                 nthreads=p.nthreads, npatterns=p.npatterns, 
                 mask_flag=p.maskflag, crop_flag=p.cropflag, 
                 nxcrop=p.nxcrop, nycrop=p.nycrop,
                 dp_shift_flag=p.dp_shift_flag, 
                 shiftx=p.shiftx, shifty=p.shifty,
                 maskname=p.path_to_string(p.maskname), rebin=p.rebin, nstart=p.nstart,
                 outputdp=p.outputdp, corrtype=p.corrtype)
    #
    # Check the rebin and crop values
    #
    if (corr.crop_flag==True) and (corr.rebin_flag==1) and ((corr.nxcrop%corr.rebin!=0) or (corr.nycrop%corr.rebin!=0)):
        print("Warning: Crop and rebin flags set, but rebin factor does not divide the nxcrop and/or nycrop values.")
        print("This may cause a sub-bin sized error in the diffraction pattern center. Maybe an issue if you have sharp Bragg peaks.")
        print("Consider values that do divide: e.g. nxcrop=1024; rebin= 2 or 4") 

    #
    # load mask and compute mask correlation 
    #
    if corr.mask_flag:
        corr.load_mask_from_file()
        corr.calculate_mask_correlation()

    #
    #  append qmax to the parameter log file
    #
    outname = p.makefname( p.outpath, p.tag, "_difftocorr_parameter_log", ".txt")
    corr.append_qmax_to_parameter_log(outname)
    p.qmax = corr.qmax

    #
    # calculate the correlation function
    #
    corrsum = np.zeros( (p.nx//2, p.nx//2, p.nth) )
    print("\nPerforming Correlations")
    print(f'Background estimate?   {corr.bg_estimate}')
    print(f'Difference correlation?  {corr.diffcorrflag}')
    print("\n")
    corrsum = corr.calculate_correlation()


    #
    # mask correction
    #   
    if corr.mask_flag:
        corrsum[:,:,:,0] = corr.mask_correction(corrsum[:,:,:,0])
        corrsum[:,:,:,1] = corr.mask_correction(corrsum[:,:,:,1])


    print("\n")
    #
    # write the correlation function to file
    #
    #outname = p.outpath+p.tag+"_correlation_sum.dbin"  #append diff or bg as appropriate
    #padfio.write_dbin( outname, corrsum ) 
    outname = p.outpath / (p.tag+"_a_correlation_sum.npy")  #append diff or bg as appropriate
    np.save( outname, corrsum[:,:,:,0] )
    print("Written correlation sum:", outname)
    if p.writeconfigs: 
        p.write_padf_config(outname, p.tag+"_a")
        p.write_mask_config(outname, p.tag+"_a")
        p.write_plot_config(outname, p.tag+"_a")

    outname = p.outpath / (p.tag+"_b_correlation_sum.npy")  #append diff or bg as appropriate
    np.save( outname, corrsum[:,:,:,1] ) 
    print("Written correlation sum:", outname)
    if p.outputsum:
        outname = p.outpath / (p.tag+"_ab_correlation_sum.npy")
        np.save( outname, corrsum[:,:,:,0]+corrsum[:,:,:,1] ) 
        print("Written correlation sum:", outname)
    #if p.writeconfigs: 
    #    p.write_padf_config(outname, p.tag+"_b")
    #    p.write_mask_config(outname, p.tag+"_b")
    #    p.write_plot_config(outname, p.tag+"_b")
    print("Correlations Done.")
