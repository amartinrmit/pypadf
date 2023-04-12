
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import subprocess
import glob

import params.paramsCORR as params
import params.paramsDIFFCORR as paramsDIFFCORR
import fxstools.diffraction as df
import fxstools.quaternions as quaternions
import fxstools.pydiffractionio as pio
import fxstools.correlation as correlation
import array
import fxstools.padfio as padfio

if __name__ == '__main__':
        #
        # set up parameter class and read all parameters from file
        #
        p = paramsDIFFCORR.paramsDIFFCORR()

        print("diffract_and_correlate.py")
        print("pypadf version ",p.version,"\n")

        #
        # Read parameters from the configuration file
        #
        p.read_config_file()
        p.qmax_calc()

        #
        # make a temporary output directory for the diffraction patterns
        #
        outpath = p.outpath
        outpathdp = outpath / "dp/"
        if not os.path.exists(p.path_to_string(outpathdp)):
            os.mkdir(p.path_to_string(outpathdp))

        #
        # Set up the diffraction class object
        #
        pdbname = p.path_to_string(p.pdbname)
        difrct = df.diffraction(pdbname, outpathdp, p.tag, p.fext, p.nx, p.wl, p.dz, p.pw, p.cenx,
                                p.ceny, p.henkeflag, np.array([p.rx,p.ry,p.rz]),p.rtheta,p.rotflag )
        p.samplepath = outpathdp

        difrct.pdb.maxdims()
        mins =  (difrct.pdb.xmin,difrct.pdb.ymin,difrct.pdb.zmin)

        if(p.alen<difrct.pdb.xlen)or(p.blen<difrct.pdb.ylen)or(p.blen<difrct.pdb.zlen):
            uc = df.unitcell(a=difrct.pdb.xlen,b=difrct.pdb.ylen,c=difrct.pdb.zlen)
        else:
            uc = df.unitcell(a=p.alen,b=p.blen,c=p.clen)

        #save the pdbname
        pdbshiftname = pdbname[:-4]+"_shifted.pdb" 

        #
        # Set up an instance of the correlation class
        #

        # make the filelist
        #
        #(read in a file list or make a list using glob)
        p.samplepath = p.path_to_string(outpathdp / (difrct.tag+"*"+difrct.fext))
        p.load_flist_from_samplepath()

        #print("DEBUG rebin:,", p.rebin)
        if p.ny<0: p.ny=p.nx   #CHECK THAT THIS CATCHES ANY ERRORS IN THE INPUT PARAMETERS?!? MAYBE DON'T RUN IF NX IS NOT SET
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


        # CHECK ALL PARAMETER NAMES GOING INTO THE DIFFRACTION AND CORRRELATION CLASSES!!!!!
        # CHECK ANY PATHS USED IN THIS MAIN SCRIPT!!!!

        #
        # calculate the 2D diffraction patterns and correlate them in chunks
        #
        ncycles = int(p.npatterns/p.chunksize)
        print("ncycles :", ncycles)
        savename1 = p.path_to_string(outpath / "SomeInitialValue1")
        savename2 = p.path_to_string(outpath / "SomeInitialValue2")

        for i in np.arange(ncycles):
            print( "\nCycle ",str(i+1)," / ",str(ncycles))
            startcycle = time.perf_counter()
            #
            # remove the chunked diffraction patterns   
            #
            #subprocess.run('rm '+p.samplepath) #might be worth making a scratch directory to ensure the wrong thing is not deleted
            
            #generate the diffraction patterns in this chunk
            dptimes = np.zeros(p.chunksize)
            for j in np.arange(p.chunksize):
                    #print("In 2nd loop", i, j)
                    if p.rotflag: difrct.axis, difrct.theta = quaternions.random_vec_angle() 

                    # shift the pdb coordinates using a unit cell
                    difrct.pdbname = pdbname
                    difrct.load_pdb()
                    if p.periodicshift: difrct.circ_shift_coordinates( uc, mins,True,p.rmax)
                    # write a new pdb file with the shifted coordinates
                    #outname = p.pdbname[:-4]+"_shifted.pdb"
                    difrct.pdbname = pdbshiftname
                    difrct.pdb.write_pdb(difrct.pdbname)
                    # load the pdb file, to ensure new atom infomration is appropriately sorted
                    difrct.load_pdb()

                    start = time.perf_counter()
                    difrct.diffraction2D()
                    end = time.perf_counter()
                    #print( "Time to calculate 2D diffraction pattern:", end-start, "seconds")
                    dptimes[j] = end-start

                    #
                    # output 2D diffraction pattern to file
                    #
                    fname = p.path_to_string( outpathdp / (difrct.tag+"_"+str(j)+difrct.fext))
                    print("DEBUG dc", fname)
                    pio.write_image( fname, difrct.dp2d )

                    #
                    # Calculate 1D pattern
                    #
                    difrct.diffraction1D()

            print( "Average time to calculate 2D diffraction patterns in this chunk:", np.average(dptimes), "seconds")

            # reload the filelist
            p.load_flist_from_samplepath()
            corr.flist = p.flist

            #
            # calculate the correlation function of the chunk
            #
            corrchunk = np.zeros( (p.nx//2, p.nx//2, p.nth) )
            print("Performing Correlations")
            corrchunk = corr.calculate_correlation()

            #
            # sum the correlations of the chunks    
            #
            if i==0:
                corrsum = corrchunk.copy()
            else:
                corrsum += corrchunk

            #
            # periodically save the correlation function
            #
            m = i%int(p.writefreq/p.chunksize)
            if m==0:
               
                #remove the previous files
                if os.path.isfile(savename1):subprocess.run( ["rm",savename1])
                if os.path.isfile(savename2):subprocess.run( ["rm",savename2])

                #save the current corrrelation file
                outname = p.path_to_string(outpath / (p.tag+"_a_npat"+str((i+1)*p.chunksize)+"_correlation_sum.npy"))  #append diff or bg as appropriate
                savename1 = outname
                np.save( outname, corrsum[:,:,:,0] ) 
                print("Written correlation sum:", outname)
                outname = p.path_to_string(outpath / (p.tag+"_b_npat"+str((i+1)*p.chunksize)+"_correlation_sum.npy"))  #append diff or bg as appropriate
                savename2 = outname
                np.save( outname, corrsum[:,:,:,1] ) 
                print("Written correlation sum:", outname)

            # cycle timing
            print( "Cycle took  :", time.perf_counter()-startcycle, " seconds")


        #
        # write the correlation function to file
        #
        #outname = p.outpath+p.tag+"_correlation_sum.dbin"  #append diff or bg as appropriate
        #padfio.write_dbin( outname, corrsum ) 
        outname = p.path_to_string( outpath / (p.tag+"_a_correlation_sum.npy"))  #append diff or bg as appropriate
        np.save( outname, corrsum[:,:,:,0] ) 
        print("Written correlation sum:", outname)
        outname = p.path_to_string( outpath / (p.tag+"_b_correlation_sum.npy"))  #append diff or bg as appropriate
        np.save( outname, corrsum[:,:,:,1] ) 
        print("Written correlation sum:", outname)
        print("Correlations Done.")



        #
        # Plot the output to screen
        #

        # q data
        nq = p.nx
        q = np.arange( nq) * 0.1 * 2 * np.pi
        invAngstromSymbol = r'$\mathrm{\AA}^{-1}$'


        gamma = 0.3
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





