"""
    Python tools to determine the shift of the centre of the diffraction pattern 
    relative to the first pattern in a series

    The relative shift is determined using a cross-correlation of each diffraction pattern 
    with a reference diffraction pattern
    
    Writes out rebinned recentred data

    Author: Andrew Martin   andrew.martin@rmit.edu.au
    2026
"""



import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import multiprocessing as mp
import glob
import os
import pathlib

def tifread(filename):
    """Reads a floating point tifffile use numpy.asarray
    """
    data = np.asarray(Image.open(filename))
    # Checks if the array is 3D then flattens it
    if len(np.shape(data)) > 2:
        dig_data = np.ones((data.shape[0], data.shape[1]))
        xlen = len(data[:, 0, 0])
        ylen = len(data[0, :, 0])
        for x in range(0, xlen):
            for y in range(0, ylen):
                dig_data[x][y] = data[x][y][0]
        return dig_data
    else:
        return data


def rebin_pattern(imagein, nbin):
        """Rebins a pattern
        """
        image = np.copy(imagein)

        imagesum = image * 0.0

        for i in np.arange(nbin) - nbin // 2:
            for j in np.arange(nbin) - nbin // 2:
                imagesum += np.roll(np.roll(image, i, 0), j, 1)
        imagesum *= 1.0 / float(nbin * nbin)
        out = imagesum[::nbin, ::nbin]
        return out


def recentre_flist( fnames, frefimage, rebinflag, rebin, outputshifted, outpath2, ithread, return_dict):
    """ Recentres all images in a file list relative to a reference image
    """
    
    if rebinflag:
        imagesum = rebin_pattern(frefimage,rebin)*0.0
    else:    
        imagesum = np.copy(frefimage)*0.0
        
    s = frefimage.shape
    for i, fname in enumerate(fnames[0:]):
        #if i>10: break
        #image  = tifread( fname ).copy()
        #image = np.load(fname)
        
        root, ext = os.path.splitext(fname)
        if (ext == '.tiff') or (ext=='.tif'):
            image = tifread( fname )
        elif ext=='.npy':
            image = np.load( fname )

        
        fimage = np.fft.fft2( image )
        cc = np.real(np.fft.ifft2( frefimage.conjugate()*fimage))
        imax = np.array(np.where( np.roll(np.roll(cc,s[0]//2,0),s[1]//2,1) == np.max(cc) )) - s[0]/2
        

        #print( i, imax )
        shifted = np.roll(np.roll( image, -int(imax[0]), 0), -int(imax[1]), 1 )
        

        outfname = os.path.basename(fname)[:-5]+".npy"
        outfile = str( outpath2 / outfname )
        if rebinflag:
            rebinned = rebin_pattern(shifted,rebin)
            np.save(outfile,rebinned)   
            imagesum += rebinned
        elif outputshifted:
            np.save(outfile,shifted)   
            imagesum += shifted
        else:
            imagesum += shifted
        #if i>300: break

    imagesum *= 1.0/len(fnames)
    return_dict[ithread]  = np.real(imagesum) 


class recenter_diffraction_patterns:
    """
    Class to multiprocess the recentring of diffraction patterns
    """

    def __init__( self,  tag, dpath="E:\\padf\\data\\", outpath="E:\\padf\\results\\", recenpath="E:\\padf\\recentred\\", rebin=4,\
                  rebinflag=True, outputshifted=True,nthreads=4, makedirs=True, fext='.tiff'):

        self.tag = tag #"tiff_scan1_c60_50GPa_Cut1_ta-C_Near_Diamond_Band_100by100pts_3by3nm_Diffraction_300mm_Alpha(3)_Spot(5)"
        self.dpath = pathlib.Path( dpath )
        self.outpath = pathlib.Path(outpath)
        if (makedirs==True) and not os.path.isdir(self.outpath.resolve() ):
            os.makedirs(self.outpath.resolve() )

        tmp = "*"+fext
        search_string = str((self.dpath / tmp).resolve())
        print( type(search_string) )
        self.fnames = glob.glob(search_string)
        print(" Number of filenames found:", len(self.fnames))
        print(search_string)

        self.rpath = self.outpath / "registration"
        if (makedirs==True) and not os.path.isdir(self.rpath.resolve()):
            os.makedirs(self.rpath.resolve() )
        tmp =  self.tag+"_registrationlist.txt"
        self.outname = str((self.rpath / tmp).resolve())
        tmp = self.tag+"imav.npy"
        self.outsumname = str((self.rpath.resolve() / tmp).resolve())

        self.outpath2 = pathlib.Path(recenpath) 
        if (makedirs==True) and not os.path.isdir(self.outpath2.resolve()):
            os.makedirs(self.outpath2.resolve() )
        
        self.rebin = rebin
        self.rebinflag = rebinflag
        self.outputshifted = outputshifted

        self.nthreads = nthreads


        

    def recentre_patterns(self, yin, xin, npatterns):
        """ recentre patterns using multithreading
        """


        if npatterns<1:
            npatterns=len(self.fnames)

        self.refname = self.fnames[0]

        root, ext = os.path.splitext(self.fnames[0])
        if (ext == '.tiff') or (ext=='.tif'):
            refimage = tifread( self.refname )
        elif ext=='.npy':
            refimage = np.load( self.refname )
        else:
            print("Unsupported file format for recentering code; Supported formats .tiff, .npy")

        frefimage = np.fft.fft2( refimage )
        s = refimage.shape

        yshift, xshift = s[0]//2-yin, s[1]//2-xin
        refimage = np.roll(np.roll( refimage, xshift,0), yshift, 1)
        frefimage = np.fft.fft2( refimage )


        manager = mp.Manager()            
        return_dict = manager.dict()
        processes = []
        for i in range(self.nthreads):
            fnames_i = self.fnames[i:npatterns:self.nthreads]

            print("Main threading loop", i)
            p = mp.Process(target=recentre_flist, args=(fnames_i,frefimage,self.rebinflag,self.rebin, self.outputshifted,self.outpath2,i,return_dict))
            p.start()
            processes.append(p)

        print("joining")
        for p in processes:
            p.join()

        
        if self.rebinflag:
            self.imagesum = np.real(rebin_pattern(refimage,self.rebin))*0.0
        else:
            self.imagesum = np.real(refimage)*0.0
        for j in np.arange(self.nthreads):
            self.imagesum += return_dict[j] 

        np.save( self.outsumname, self.imagesum)


