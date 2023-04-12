"""
blqq.py

Store B_l(q,q') arrays and compute required number of spherical Bessel zeros

Classes:
    Blqq
    blqqarray
"""

import numpy as np
import fxstools.vol as vol
from scipy.special import spherical_jn
import os
import sys

class Blqq:
    """
    Stores a single B_l(q,q') matrix
    
    Attributes
    ----------
    l : int
        spherical harmonic order of the B_l(q,q') matrix
    
    nq : int
        number of radial (q) bins

    data : numpy array (float)
        B_l(q,q') data of size nq x nq
    """
    def __init__(self,l, nq):
        """Construct the Blqq class
        """
        self.l = l
        self.nq = nq
        self.data = np.zeros( (nq,nq))



class blqqarray:
    """
    Store all the Blqq arrays up to a maximum l value
 
    Attributes
    ----------
    nl : int
        maximum spherical harmonic order

    nq : int
        number of q bins

    qmax : float
        q-value of largest q bin

    rmax : float
        r-value of largest real-space radial (r) bin

    tablenzero : int
        number of n values in spherical Bessel zero lookup table
    
    tablelmax : int
        number of l values in spherical Bessel zero lookup table
    
    mult : float
        multiply max zero condition (to increase number of zeros used)
    """

    def __init__(self, nl, nq=-1, qmax=-1, rmax=-1,tablenzero=1000,tablelmax=100,mult=1):
        """
        Constructs a class to store all Blqq arrays
 
        Parameters
        ----------
        nl : int
            maximum spherical harmonic order

        nq : int
            number of q bins

        qmax : float
            q-value of largest q bin

        rmax : float
            r-value of largest real-space radial (r) bin

        tablenzero : int
            number of n values in spherical Bessel zero lookup table
        
        tablelmax : int
            number of l values in spherical Bessel zero lookup table
        
        mult : float
            multiply max zero condition (to increase number of zeros used)
        """

        
        self.nl = nl
        self.nq = nq
        self.qmax = qmax
        self.rmax = rmax
        self.l = []
        self.mult = mult

        #
        # parameters to read a lookup table of spherical bessel zeros
        #
        self.tablenzero = tablenzero
        self.tablelmax = tablelmax
        self.read_sphB_zeros()

        #
        # set up the sphB matrices
        #
        if (qmax>0) and (rmax>0):
            self.set_up_blqq_array_sphB()
        elif (nq>0):
            self.set_up_blqq_array()
        else:
            print("Error - blqq array cannot be initialised. Check qmax, rmax or nq values.")
            exit()


    def set_up_blqq_array(self):
        """
        Create an attribute (l) that lists all of the blqq objects 
        with regular q sampling
        """
        self.l = []
        for i in np.arange(self.nl):
            self.l.append( Blqq(i,self.nq) )


    def set_up_blqq_array_sphB(self):
        """        
        Create an attribute (l) that lists all of the blqq objects 
        with spherical Bessel zero sampling
        """

        self.l = []
        for i in np.arange(self.nl):
            nq = self.sphB_samp_nmax( i, self.rmax, self.qmax )
            self.l.append( Blqq(i,nq) )


    def write_blqq_array_as_dbin(self,path,tag):
        """
        Write all the blqq arrays out to dbin files

        Parameters
        ----------
        path : str
            output path where files will be written

        tag : str
            prefix for the start of each file name
        """

        for i in np.arange(self.nl):
            if type(path)=='str':
                outname = path+tag+"_"+str(i)+".dbin"
            else:
                outname = path / (tag+"_"+str(i)+".dbin")
                outname = str(outname.resolve())
            padfio.write_dbin( outname, self.l[i].data )

    def write_blqq_array_as_npy(self,path,tag):
        """
        Write all the blqq arrays out to npy files

        Parameters
        ----------
        path : str
            output path where files will be written

        tag : str
            prefix for the start of each file name
        """
        for i in np.arange(self.nl):
            if type(path)=='str':
                outname = path+tag+"_"+str(i)+".npy"
            else:
                outname = path / (tag+"_"+str(i)+".npy")
                outname = str(outname.resolve())
            np.save( outname, self.l[i].data )


    def sphB_samp_nmax( self, l, rmax, qmax):

        """
        Compute the number of spherical Bessel zeros up to 2pi*rmax*qmax

        Parameters
        ----------
        l : int
            spherical harmonic order to compute

        rmax : float
            maximum real-space radial value

        qmax : float
            maximum q-space radial value

        Returns
        -------
        out : int
            number of spherical Bessel zeros
        """
        
        qlim = 2*np.pi*qmax*rmax*self.mult
        out= -1
        for i in np.arange(self.tablenzero):
            qln = self.jnuzeros[ l, i ]
            #print("l i qln qlim", l, i, qln, qlim)
            if qln>qlim :
                out = i-1
                break
        if out<0:
            out=0

        return out       

    def read_sphB_zeros(self):
        """
        Read a table of spherical Bessel zero values from a file
        """
        tablename = os.path.join( sys.path[0], "sphbzeros_lmax"+str(self.tablelmax)+"_nt"+str(self.tablenzero)+".npy" )
        self.jnuzeros = np.load(tablename)


    #
    # resample from zerom from one order to zeros of another
    #
    def resampling_matrix(self,l, l2, nmax, nmax2, qmax, rmax ):
       """
        Compute the resampling matrix to interpolate from zeros spherical harmonic 
        of order l to zeros of order l2

        Parameters
        ----------
        l : int
            input spherical harmonic order

        l2 : int
            target spherical harmonic order

        nmax : int
            number of radial bins

        nmax2 : int
            number of radial bins

        qmax : float
            q value of largest radial bin

        rmax : float
            r value of largest real-space radial bin.

        Returns
        -------
        mat : numpy array (float)
            resampling matrix
       """
        
       if (qmax<0) or (rmax<0):
           print("Error - resampling matrix calc with negative rmax/qmax vals")
           exit()

       #print( "resampling matrix, rmax, qmax", rmax, qmax )

       qmax2p = 2*np.pi*qmax

       mat = np.zeros( (nmax2, nmax) )

       q2arr = np.copy(self.jnuzeros[l2, 1:nmax2+1 ]/rmax)

       qarr  = np.copy(self.jnuzeros[l, 1:nmax+1])
       jl2arr = spherical_jn( l+1, qarr )
       qarr *= 1.0/rmax


       r = np.copy(self.jnuzeros[l, 1:nmax+1] / qmax2p)
       jl1 = jl2arr


       for i in np.arange(nmax2):
           q2 = q2arr[i]

           for j in np.arange(nmax):
               q = qarr[j]
               jl2 = jl2arr[j]

               factor2 = np.sqrt(2*np.pi)/ ( jl2*jl2*(rmax**3))
               sum = 0.0
               for k in np.arange(nmax):
                   
                   factor = np.sqrt(2*np.pi)/ ( jl1[k]*jl1[k]*(qmax2p**3))                   
                   sb1 = spherical_jn( l, q*r[k])
                   sb2 = spherical_jn( l, q2*r[k])
                   sum += sb1*sb2*factor*factor2

               mat[i,j] = sum

       return mat



               
