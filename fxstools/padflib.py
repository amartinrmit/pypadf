"""
padflib.py

Tools to calculate the PADF from the correlation function
and related quantities (e.g. B_l(q,q') matrices)
Based on the discrete spherical Bessel transform
sampled at the spherical Bessel zeros.
(include citation)

Classes:
    padfcl - padf/B_l(q,q') calculation
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.special import spherical_jn, legendre
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
import scipy as sp
import time

import fxstools.vol as vol
import fxstools.thmask as thmask
import fxstools.blqq as blqq


#
# TO DO:
#   -  check odd/even parts
#   -  check all the matrix multiplcations
#   -  add timings ans outputs
#   -  outputs after each step
#
#   - write a test script; uses new parameter class and new padf class
#     imports the old results from c-code and checks all the outputs
#
#   - after that all works; continue onto the correlation parts of the code

class padfcl:
    """
    PADF and B_l(q,q') calculation tools

    Parameters
    ----------
    nl : int
        maximum order of spherical harmonics (l)

    nlmin : int
        minimum spherical harmonic order to use

    nr : int
        number of radial bins in real space

    nth : int
        number of angular bins

    nq : int
        number of radial bins in q-space
    
    qmin : float
        minimum value of q in the correlation file. Units: inverse metres
        (NOTE q IS DEFINED WITHOUT A FACTOR OF 2PI CONTRARY TO USUAL X-RAY CONVENTIONS)

    qmax : float
        maximum value of q in the correlation file. Units: inverse metres 
        (NOTE q IS DEFINED WITHOUT A FACTOR OF 2PI CONTRARY TO USUAL X-RAY CONVENTIONS)
    
    rmax : float
        maximum value of r in the correlation file. Units: metres

    wl : float
        wavelength of the incident beam. Units: metres

    corrvol : array-like (floats)
        floating point valued numpy array that stores the q-space correlation volume

    units_flag : bool
        rescale the padf volume to physical units (not implemented)

    tablenzero : int
        maximum n value in spherical harmonic look-up table

    tablelmax : int
        maxmimy l value in spherical harmonic look-up table

    method : str
        matrix inversion method for obtaining B_l matrices from correlation volume
        'svd' - singular value decomposition
        'legendre' - use orthogonality properties of Legendre polynomials
    
    mult : int 
        multiplication factor to increase the number of spherical Bessel zeros
        included in the calcualtion (default value = 1)

    legendre_norm : bool
        Normalise the Legendre polynomials. 
        True is numerically accurate.
        False is consistent with earlier padf versions.

    padf : vol object
        the final padf volume

    blqq : arrary of blqq objects
        the blqq arrays for all orders of l
    """

    
    def __init__(self, nl=2, nlmin=0, nr=100, nth=100, nq=100, qmin=0.0, qmax=0.0, rmax=0.0, corrvol=None,\
                    wl=0.0, units_flag=0,\
                    tablenzero=1000, tablelmax=100, method='svd',mult=1, legendre_norm=True):
        """
        Constructs an instance of the padfcl class
        Initialises padf and blqq objects
        """ 

        self.nl = nl
        self.nlmin = nlmin
        self.nr = nr
        self.nth = nth
        self.nq = nq
        self.qmin = qmin
        self.qmax = qmax
        self.rmax = rmax
        self.corrvol = corrvol
        self.wl = wl
        self.units_flag = units_flag # not implemented
        self.method=method 
        self.mult = mult
        self.legendre_norm = legendre_norm

        #
        # parameters to read a lookup table of spherical bessel zeros
        #
        self.tablenzero = tablenzero
        self.tablelmax = tablelmax
        

        #
        # creates and empty volume object
        #
        dimnames = [0,0,0] #"q","q2","theta"]
        dimlen = [self.nq,self.nq,self.nth]
        dmin = [0,0,0]
        dmax = [self.qmax,self.qmax,2.0*np.pi] 

        self.padf = vol.Vol2(  dimnames=dimnames, dimlen=dimlen, dmin=dmin, dmax=dmax )

        #print("DEBUG padflibdev qmax, rmax:", self.qmax, self.rmax)
        #print("DEBUG PADFLIBDEV tablenzero tablelmax", self.tablenzero, self.tablelmax)
        self.blqq = blqq.blqqarray(self.nl, qmax=self.qmax, rmax=self.rmax,
                                  tablenzero=self.tablenzero,tablelmax=self.tablelmax,mult=self.mult) 



    #
    #  The main calculation of the padf from the correlation volume
    #
    def padf_calc( self ):
       """
        Calculate the padf volume from the correlation function.
        All key parameters must be set when padfcl in constructed.
       """
       #print("DEBUG padflibdev.py "+str(self.legendre_norm))
       print("\nPADF calculation started...")
       start = time.perf_counter()
       #
       # STEP 1: calculate the B_l(q,q') matrices sampled at the spherical Bessel zeros
       #
       #self.Blqq_calc_fast(method='svd')
       #print( "DEBUG padf_calc self.method=", self.method)
       #self.Blqq_calc_fast(method=self.method,outputfmat=True)
       self.Blqq_calc_fast_interp(method=self.method,order=1)
       tpoint2 = time.perf_counter()
       print("Total time to calculate step 1 (correlation->blrr):", tpoint2-start," seconds")

       #
       # STEP 2: transform the B_l(q,q') into B_l(r,r') terms
       #
       self.blrr = blqq.blqqarray(self.nl, nq=self.nr ) 
       self.Bla_qr_transform()
       
       tpoint3 = time.perf_counter()
       print("Total time to calculate step 2 (blqq->blrr):", tpoint3-tpoint2," seconds")
       #
       # STEP 3: transform B_l(r,r') into the PADF
       #
       self.padf = self.Blrr_to_padf( self.blrr.l, (self.nr,self.nr,self.nth) )
       end = time.perf_counter()
       print("Total time to calculate step 3 (blrr->padf):", end-tpoint3," seconds")
       print("Total time to calculate PADF:", end-start," seconds")

    def normalise_padf_with_density( self, density):
       """make the padf dimenionsionless by dividing by the density in m^-3
       """ 
       self.padf *= 1.0/(density**4) #


    #
    # test the blqq_to_corr function 
    #
    def test_blqq_to_corr( self ):
       """
        Test the correlation->Blqq operation.
        Calculates the Blqq matrices from a correlation function
        then recalculates the correlation function from the blqq arrays.
        Print a sum square difference of the initial and final correlation volumes.
       """

       #
       # STEP 1: calculate the B_l(q,q') matrices sampled at the spherical Bessel zeros
       #
       self.Blqq_calc_fast()
       #self.Blqq_calc_fast_interp()

       #self.blqq = self.Blqq_Ewald_filter( self.blqq, self.corrvol.shape )
       
       # Now test going back the other way
       corr2 = self.Blqq_to_corr_fast(self.blqq)

       print("corr vol difference:", np.sqrt(np.sum( (corr2-self.corrvol)**2)), 
               np.sqrt(np.sum( (corr2)**2)), 
               np.sqrt(np.sum( (self.corrvol)**2)) )
       self.corrvol = corr2

    #
    # Generates self.blqq matrices from the correlation function
    # Each matrix is sampled at the zeros the corresponding spherical bessel function
    #
    def Blqq_calc_fast(self, method='svd', outputfmat=False):
        """
        Compute the B_l(q,q') arrays from the correlation volume.
        The q-sampling for the l=0 spherical Bessel zeros are used.
        This 'fast' because the same q-sampling is used for each order l.
        
        Parameters
        ---------
        method : str
            matrix inversion method for obtaining B_l matrices from correlation volume
            'svd' - singular value decomposition
            'legendre' - use orthogonality properties of Legendre polynomials
    
        outputfmat : bool
            save the F-matrix singular values to file

        Returns
        -------
            N/A
            Result stored in self.blqq and self.blqqtmp attributes
        

        """
        self.blqqtmp = blqq.blqqarray(self.nl, nq=self.blqq.l[0].nq, 
                                      tablenzero=self.tablenzero,tablelmax=self.tablelmax,mult=self.mult) 
        
        # define a list of matrices to store the Blqq matrices
        rsinvlist = []
        
        ###plt.figure()
        for l in np.arange(self.nl):
            mat = self.blqq.resampling_matrix(l, 0, self.blqq.l[l].nq, self.blqq.l[0].nq,\
                                                 self.qmax, self.rmax )
            
            u, s, vh = np.linalg.svd( mat, full_matrices=False )
            #print( "resamp smax smin", s[0], s[-1], "; dims u s vh", u.shape, s.shape, vh.shape )
            igood = np.where(s > 0.5)
            sinv = s*0.0
            sinv[igood] = 1.0/s[igood]
            inv = np.dot(np.dot(vh.transpose(), np.diag(sinv)), u.transpose())
            #print("resamp inv shape", inv.shape, np.max(inv), np.min(inv) )
#            inv = sp.linalg.pinv(mat)    # I MAY NOT HAVE TO INVER THIS, BUT JUST SOLVE STUFF
            rsinvlist.append(inv)

            ###plt.plot(s)
        ###plt.draw()
        ###plt.show()
            
        nq = self.blqq.l[0].nq
        #print("DEBUG BLQQ_CALC_FAST nq", nq)
        self.fsing = []
        for iq in np.arange(nq):
           for iq2 in np.arange(nq):

              qln1 = self.blqq.jnuzeros[0,iq+1] / (2*np.pi*self.rmax*self.mult)
              qln2 = self.blqq.jnuzeros[0,iq2+1] / (2*np.pi*self.rmax*self.mult)
              
              ic = np.round( self.nq*(qln1/self.qmax) ).astype(np.int)
              jc = np.round( self.nq*(qln2/self.qmax) ).astype(np.int)
              if (ic>=nq) or (jc>=nq): continue 
              
              # CALCULATE AND INVERT FMAT; MAY NOT HAVE TO INVERT JUST SOLVE
              if method=='svd':
                 mat = self.fmat( qln1, qln2 )
                 zmat = self.fmat_zsamp( qln1, qln2)
                 if outputfmat:
                     #testmat = mat
                     testmat = zmat #self.fmat_zsamp( qln1, qln2 )
                     normmat = np.dot( testmat, testmat.transpose() )

                     u, s, vh = np.linalg.svd( testmat )
                     self.fsing.append( s )
                     if iq==0 and iq2==0: np.save( "testmat.npy", testmat )

                     #print("DEBUG Pn orthog:", np.sum(testmat[0,:]*testmat[1,:] ))

                 #matinv = sp.linalg.pinv2( mat, cond=0.5)   
                 #tmp = np.dot( matinv.transpose(), self.corrvol[ic,jc,:])
                 
                 
                 z = np.arange(self.nth)*2/self.nth - 1
                 acs = (self.nth/2)*np.arccos(z)/np.pi
                 acs2 = (self.nth/2)*np.arccos(z[::-1])/np.pi + self.nth/2
                 zorder = 1
                 zcorr = map_coordinates(self.corrvol[ic,jc,:], [acs], order=zorder) 
                 zcorr +=  map_coordinates(self.corrvol[ic,jc,:], [acs2], order=zorder) 

                 zmatinv = sp.linalg.pinv( zmat, cond=0.5 )
                 tmp = np.dot( zmatinv.transpose(), zcorr )
                 
              elif method=='legendre':
                 tmp = self.Blqq_legendre( qln1, qln2, self.corrvol[ic,jc,:], True )

              else:            
                 print("invalid method option for calculating Blqq. 'svd' or 'legendre'. Exiting.")
                 exit()
 
              for k in np.arange(self.nl//2):
                 self.blqqtmp.l[2*k].data[iq,iq2] = tmp[k]
                 
              for k in np.arange(self.nl//2-1):
                 self.blqqtmp.l[2*k+1].data[iq,iq2] = 0.0

        # RESAMPLE USING THE RESAMPLING MATRICES...
        for k in np.arange(self.nl//2):
           tmp = np.dot(rsinvlist[2*k],self.blqqtmp.l[2*k].data)
           self.blqq.l[2*k].data = np.dot( tmp, rsinvlist[2*k].transpose())

        if outputfmat:
             #print("DEBUG output fmat singular values")
             self.fsing = np.array(self.fsing)
             np.save( "fsing.npy", self.fsing )   
 
    #
    # Generates correlation volume from the blqq matrices 
    # Each matrix is sampled at the zeros the corresponding spherical bessel function
    #
    def Blqq_to_corr_fast(self,blqq_in,interpolate=False, order=3):
        """
        Compute the correlation volume from the B_l(q,q') arrays
        using the spherical Bessel zero sampling.
        
        Parameters
        ---------
        blqq_in : blqq object
            list of blqq array objects

        interpolate : bool
            use spline interpolation to rebin q-sampling
            If False, use nearest neighbour sampling
    
        order : int
            order for spline interpolation
            used if interpolate==True

        Returns
        -------
        corrvol : array like (floats)
            q-space correlation volume

        """

        corrvol = np.zeros(( self.nq, self.nq, self.nth))

        self.blqqtmp = blqq.blqqarray(self.nl, nq=blqq_in.l[0].nq, 
                                      tablenzero=self.tablenzero,tablelmax=self.tablelmax) 
        
        # define a list of matrices to store the Blqq matrices
        rsinvlist = []
        rslist = []

        for l in np.arange(self.nl):
            mat = self.blqq.resampling_matrix(l, 0, blqq_in.l[l].nq, blqq_in.l[0].nq,\
                                                 self.qmax, self.rmax )
             
            rslist.append(mat)
            
            
        for k in np.arange(self.nl//2):
           tmp = np.dot(rslist[2*k],blqq_in.l[2*k].data)
           self.blqqtmp.l[2*k].data =  np.dot( tmp, rslist[2*k].transpose())
        
        zlist = []
        nq = blqq_in.l[0].nq
        #print("DEBUG BLQQ_TO_CORR_FAST nq", nq)
        corrvoltmp = np.zeros( (nq,nq,self.nth))
        for iq in np.arange(nq):
           for iq2 in np.arange(nq):

             qln1 = blqq_in.jnuzeros[0,iq+1] / (2*np.pi*self.rmax*self.mult)
             qln2 = blqq_in.jnuzeros[0,iq2+1] / (2*np.pi*self.rmax*self.mult)
              
             # CALCULATE AND INVERT FMAT; MAY NOT HAVE TO INVERT JUST SOLVE
             mat = self.fmat( qln1, qln2 )
             #matinv = sp.linalg.pinv2( mat, cond=0.5)
              

             if interpolate==True:
                  for k in np.arange(self.nl//2):
                     corrvoltmp[iq,iq2,:] += mat[k]*self.blqqtmp.l[2*k].data[iq, iq2]
             
             else:
                  # FIND THE NEAREST Q INDEX IN THE CORRELATION DATA
                  ic = np.round( self.nq*(qln1/self.qmax) ).astype(np.int)
                  jc = np.round( self.nq*(qln2/self.qmax) ).astype(np.int)
                  
                  if (ic>=nq) or (jc>=nq): continue 
                 
                  for k in np.arange(self.nl//2):
                     corrvol[ic,jc,:] += mat[k]*self.blqqtmp.l[2*k].data[iq, iq2]
         
        if interpolate:
            for iq in np.arange(self.nq):
                for iq2 in np.arange(self.nq):
                      fact = (2*np.pi*self.qmax*self.rmax*self.mult)/self.nq
                      zdist = np.sqrt( (iq*fact - blqq_in.jnuzeros[0,:])**2 )
                      iz = np.where(zdist == np.min(zdist) )
                      zdist2 = np.sqrt( (iq2*fact - blqq_in.jnuzeros[0,:])**2 )
                      iz2 = np.where(zdist2 == np.min(zdist2) )
                      #print( "interpolate test ", iz[0], iq*fact, iz2[0], iq2*fact )
                      zlist.append( np.array( [iz[0][0], iz2[0][0]] ) )                            
            zlist = np.array(zlist) 
            for ith in np.arange( self.nth ):
                 #print("\r test ith", ith)
                 corrvol[:,:,ith] = map_coordinates( corrvoltmp[:,:,ith], zlist.transpose(), order=order ).reshape( (self.nq,self.nq) )

            #print("DEBUG corrvol max", np.min(corrvol), np.max(corrvol), np.min(corrvoltmp), np.max(corrvoltmp), np.min(zlist), np.max(zlist) )        
        return corrvol


    def Blqq_ewald_filter(self,blqq_in,thmin=0,qwid=0):
        """
        Filters the Blqq matrices using the Blqq->correlation transformation.
        A step in modifying an ideal padf volume to compare to an experimental padf.
        

        Parameters
        ---------
        blqq_in : blqq object
            list of blqq array objects
    
        thmin : float
            minimum value of theta in degrees

        qwid : float
            pixel width of a radial-q Gaussian filter
            if qwid==0, no Gaussian filter is applied

        Returns
        -------
        blqqout : array of blqq objects
            filtered blqq arrays    
        """

        corrvol = np.zeros(( self.nq, self.nq, self.nth))

        blqqtmp = blqq.blqqarray(self.nl, nq=blqq_in.l[0].nq, 
                                      tablenzero=self.tablenzero,tablelmax=self.tablelmax) 
        blqqtmp2 = blqq.blqqarray(self.nl, nq=blqq_in.l[0].nq, 
                                      tablenzero=self.tablenzero,tablelmax=self.tablelmax) 
        
        blqqout = blqq.blqqarray(self.nl, qmax=self.qmax, rmax=self.rmax,       
                                  tablenzero=self.tablenzero,tablelmax=self.tablelmax,mult=self.mult) 
        # define a list of matrices to store the Blqq matrices
        rsinvlist = []
        rslist = []

        for l in np.arange(self.nl):
            mat = self.blqq.resampling_matrix(l, 0, blqq_in.l[l].nq, blqq_in.l[0].nq,\
                                                 self.qmax, self.rmax )
            rslist.append(mat)
            
            u, s, vh = np.linalg.svd( mat, full_matrices=False )
            #print( "resamp smax smin", s[0], s[-1], "; dims u s vh", u.shape, s.shape, vh.shape )
            igood = np.where(s > 0.5)
            sinv = s*0.0
            sinv[igood] = 1.0/s[igood]
            inv = np.dot(np.dot(vh.transpose(), np.diag(sinv)), u.transpose())
            #print("resamp inv shape", inv.shape, np.max(inv), np.min(inv) )
#            inv = sp.linalg.pinv(mat)    # I MAY NOT HAVE TO INVER THIS, BUT JUST SOLVE STUFF
            rsinvlist.append(inv)


        for k in np.arange(self.nl//2):
           tmp = np.dot(rslist[2*k],blqq_in.l[2*k].data)
           blqqtmp.l[2*k].data =  np.dot( tmp, rslist[2*k].transpose())
        
        #print("blqqtmp (1) max", np.max(blqqtmp.l[2].data))

        nq = blqq_in.l[0].nq
        #print("DEBUG BLQQ_TO_CORR_FAST nq", nq)
        for iq in np.arange(nq):
           print("Calculating Amatrix filter iq = ",iq,"/",nq-1)
           for iq2 in np.arange(nq):

              qln1 = blqq_in.jnuzeros[0,iq+1] / (2*np.pi*self.rmax)
              qln2 = blqq_in.jnuzeros[0,iq2+1] / (2*np.pi*self.rmax)
              
              thq = self.thetaq_calc( qln1/(2*np.pi), self.wl )
              thq2 = self.thetaq_calc( qln2/(2*np.pi), self.wl )
              phi = 2 * np.pi * np.arange(self.nth)/self.nth
              arg = np.cos(thq)*np.cos(thq2) + np.sin(thq)*np.sin(thq2)*np.cos(phi)
              arg[arg>1] = 1.0
              arg[arg<-1] = -1.0
              costhmask = thmask.make_costheta_mask_v2( self.nth, -1, np.min(arg)  )
              if thmin>0: costhmask *= thmask.make_costheta_mask( self.nth, 0, thmin )

              Amatrix = thmask.costheta_mask_sphB( self.nl, costhmask )
              
              #print("Amatrix max", iq, iq2, np.max(Amatrix))

                # FIND THE NEAREST Q INDEX IN THE CORRELATION DATA
              #ic = np.round( self.nq*(qln1/self.qmax) ).astype(np.int)
              #jc = np.round( self.nq*(qln2/self.qmax) ).astype(np.int)
              #
              #if (ic>=nq) or (jc>=nq): continue 
           
              if qwid>0:
                 fact = np.exp( -0.5*(iq/qwid)**2) * np.exp(-0.5*(iq2/qwid)**2)
              else:
                 fact = 1.0
              print( iq, iq2, "fact", fact)
              for k in np.arange(self.nl//2):
                 for k2 in np.arange(self.nl//2):
                     if qwid>0:
                         fact = np.exp( -0.5*(k/qwid)**2) * np.exp(-0.5*(k2/qwid)**2)
                     else:
                         fact = 1.0

                     blqqtmp2.l[2*k].data[iq,iq2] += Amatrix[2*k,2*k2]*blqqtmp.l[2*k2].data[iq,iq2]*fact
    

        print("blqqtmp2 (2) max", np.max(blqqtmp2.l[2].data))

        # RESAMPLE USING THE RESAMPLING MATRICES...
        for k in np.arange(self.nl//2):
           tmp = np.dot(rsinvlist[2*k],blqqtmp2.l[2*k].data)
           blqqout.l[2*k].data = np.dot( tmp, rsinvlist[2*k].transpose())
         
        return blqqout

    #
    # calculate F-matrix; for conversion of correlation volume to B_l(q,q') matrices
    #
    def fmat( self, q, q2 ):
       """        
        Calculates F-matrix for a particular values of q and q'. 
        Matrix that maps B_l(q,q') into a corrrelation volume.
        Accounts for Ewald sphere curvature.
        Uses regular theta sampling.

        Parameters
        ---------
        q, q2 : float
            evaluate F-matrix at these q values. Units: inverse metres

        Returns
        -------
        mat : numpy array object (float)
            F matrix evaluated at q and q2. 
       """
       thq = self.thetaq_calc( q, self.wl )
       thq2 = self.thetaq_calc( q2, self.wl )

       phi = 2 * np.pi * np.arange(self.nth)/self.nth
       arg = np.cos(thq)*np.cos(thq2) + np.sin(thq)*np.sin(thq2)*np.cos(phi)
       arg[arg>1] = 1.0
       arg[arg<-1] = -1.0
       
       mat = np.zeros( (self.nl//2,self.nth) )
       for l in np.arange( self.nl//2 ):
          Pn = legendre( int(2*l) )
          if self.legendre_norm:
            mat[l,:] =  Pn(arg)*np.sqrt(2*(2*l)+1)
          else:
            mat[l,:] =  Pn(arg)
    
       return mat

    def fmat_zsamp( self, q, q2 ):
       """        
        Calculates F-matrix for a particular values of q and q'. 
        Matrix that maps B_l(q,q') into a corrrelation volume.
        Accounts for Ewald sphere curvature.
        Uses regular sampling of cos(theta).

        Parameters
        ---------
        q, q2 : float
            evaluate F-matrix at these q values. Units: inverse metres

        Returns
        -------
        mat : numpy array object (float)
            F matrix evaluated at q and q2. 
       """
       thq = self.thetaq_calc( q, self.wl )
       thq2 = self.thetaq_calc( q2, self.wl )
       

       #phi = 2 * np.pi * np.arange(self.nth)/self.nth
       z = np.arange( self.nth)*2/self.nth - 1
       arg = np.cos(thq)*np.cos(thq2) + np.sin(thq)*np.sin(thq2)*z
       #arg = np.cos(thq)*np.cos(thq2) + np.sin(thq)*np.sin(thq2)*np.cos(phi)
       arg[arg>1] = 1.0
       arg[arg<-1] = -1.0
       #print("DEBUG z-arg ", np.sum(np.abs(z-arg) ) )  
       mat = np.zeros( (self.nl//2,self.nth) )
       for l in np.arange( self.nl//2 ):
          Pn = legendre( int(2*l) )
          if self.legendre_norm:
            mat[l,:] =  Pn(arg)*np.sqrt(2*(2*l)+1)
          else:
            mat[l,:] =  Pn(arg)

       return mat

    def Blqq_legendre( self, q, q2, corrline, applyLimits=True):
        """       
        Computes Blqq matrix values at radial q values: q and q2
        using Lengendre function orthogonality properties
 
        Parameters
        ---------
        q, q2 : float
            radial q values

        corrline : numpy array (float)
            1D array values of correlation function at q, q2 for all theta values

        appyLimits : bool
            enforces minimum and maximum theta values

        Returns
        -------
        tmp : numpy array (float)
            1D array of B_l(q,q2') values as function of l, evaluated at q and q2. 
        """
        thq = self.thetaq_calc( q, self.wl )
        thq2 = self.thetaq_calc( q2, self.wl )
         
 
        phi = 2 * np.pi * np.arange(self.nth)/self.nth
        arg = np.cos(thq)*np.cos(thq2) + np.sin(thq)*np.sin(thq2)*np.cos(phi)
        arg[arg>1] = 1.0
        arg[arg<-1] = -1.0
        minarg, maxarg = np.min(arg), np.max(arg)  
 
        f = interp1d( arg, corrline )

        #z = 2*np.arange(self.nth)/float(self.nth) - 1
        z = (maxarg-minarg)*np.arange(self.nth)/self.nth - (maxarg-minarg)/2
        z[z<minarg] = minarg + 1e-4
        ###print( "max/min args", minarg, maxarg, np.min(z), np.max(z) )
      
        corr_z = f(z)
        if applyLimits:
          corr_z[z<minarg] = 0.0
          #corr_z[z>maxarg] = 0.0

        tmp = np.zeros( self.nl )
        for l in np.arange( self.nl//2 ):
           Pn = legendre( int(2*l) )
           if self.legendre_norm:
                factor = np.sqrt(2*(2*l)+1)
                print("test norm")
           else:
                print("test no norm")
                factor = 1
           tmp[l] =  np.sum( Pn(z)*factor * corr_z )

        return tmp
    
    #
    # Calculate the 2 x scattering angle from the magnitude of the q vector
    #
    def thetaq_calc( self, q, wl ):
        """ 
        Computes 2\theta scattering angle from q and wavelength        
       
        Parameters
        ---------
        q : float
            radial q value. Units: inverse metres

        wl : float
            wavelength. Units: metres


        Returns
        -------
        scattering angle (2 \theta)
        """
        arg = wl*q/2
        if arg>1: arg=1
        if arg<-1: arg=-1

        return  (np.pi/2)-np.arcsin(arg)



    #
    # transform blqq into real-space blrr matrices
    #
    def Bla_qr_transform(self):
        """
        Computes B_l(r,r') matrices from B_l(q,q') matrices.
        Applies the spherical Bessel transform to each q dimension.
        Uses self.blqq as input
        Writes output to self.blrr
        """
        for k in np.arange(self.nl//2):
           l = 2*k
           mat = self.dsBmatrix( l, self.blqq.l[l].nq )
           tmp = np.dot(mat,self.blqq.l[l].data)
           self.blrr.l[l].data = np.dot( tmp, mat.transpose())
          
           #matinv = self.dsBmatrix_inv(l, self.blqq.l[l].nq, 0.0)
           #tmp = np.dot(matinv.transpose(),self.blqq.l[l].data)
           #self.blrr.l[l].data = np.dot( tmp, matinv)



    #
    # Calculate the matrix for the q->r transform
    #
    def dsBmatrix( self, l, nq ):
        """        
        Computes the spherical Bessel tranform matrix for harmonic order l

        Parameters
        ---------
        l : int
            order of spherical harmonic
    
        nq : int
            number of radial q bins

        Returns
        -------
        dsbmat : numpy array (float)
            spherical Bessel transform matrix
        """
        dsbmat = np.zeros( (self.nr, nq) )
        r = np.arange(self.nr)*self.rmax/self.nr
        
        for j in np.arange(nq):
           qln = self.blqq.jnuzeros[ l, j+1 ]
           arg = qln*r/self.rmax
           jl1 = spherical_jn( l+1, qln )
           factor= np.sqrt(2*np.pi) / (jl1*jl1*self.rmax**3)
           dsbmat[:,j] = spherical_jn( l, arg )* factor
 
        return dsbmat

    #
    # r->q matrix transform
    #
    def dsBmatrix_inv(self,l, nq, qmin):       
        """        
        Computes the inverse sherical Bessel transform matrix. 

        Parameters
        ---------
        l : int
            order of spherical harmonic
    
        nq : int
            number of radial q bins

        qmin : float
            minimum q value

        Returns
        -------
        dsbmatinv : numpy array (float)
            inverse spherical Bessel transform matrix
        """
       
        dsbmatinv = np.zeros( (self.nr, nq) )
        r = np.arange(self.nr)*self.rmax/self.nr
       
        for j in np.arange(nq):
           qln = self.blqq.jnuzeros[ l, j+1 ]
           arg = qln*r/self.rmax
           jl1 = spherical_jn( l+1, qln )
           factor= np.sqrt(2*np.pi) #/ (jl1*jl1*self.rmax**3)
          
           if qln>(qmin*2*np.pi*self.rmax):
              dsbmatinv[:,j] = spherical_jn( l, arg )* np.sqrt(2/np.pi)*r*r*self.rmax*1e30/self.nr
           else:
              dsbmatinv[:,j] = 0.0

        return dsbmatinv.transpose()


    #
    # Convert the blrr matrices to the PADF (l->theta)
    #
    def Blrr_to_padf( self, blrr, padfshape ):
        """        
        Transforms B_l(r,r') matrices to padf volume
        using Legendre polynomials        

        Parameters
        ---------
        blrr : blqq object
            input blrr matrices

        padfshape : tuple (floats)
            shape of padf array

        Returns
        -------
        padfout : vol object
            padf volume
        """
        padfout = np.zeros( padfshape )
        lmax = len(blrr)
        for l in np.arange(2,lmax):

             if (l%2)!=0:
                  continue

             s2 = padfout.shape[2]
             z = np.cos( 2.0*np.pi*np.arange(s2)/float(s2) )
             Pn = legendre( int(l) )
             #print("Blrr to padf (1)"+str(l)) #DEBUG
             if self.legendre_norm:
                p = Pn(z)*np.sqrt(2*l+1)
                #print("blrr to padf; lnorm has been done")
             else:
                p = Pn(z)
                #print("lnorm has not been done")
      
             for i in np.arange(padfout.shape[0]):
                  for j in np.arange(padfout.shape[1]):
                       padfout[i,j,:] += blrr[l].data[i,j]*p[:]

        return padfout

    def project_padf_Legendre( self, padf, l ):
        """        
        Compute B_l(r,r') matrix from padf for a particular value of l
        Uses regular theta sampling.

        Parameters
        ---------
        padf : vol object
            input padf volume

        l : int
            spherical harmonic order

        Returns
        -------
        output : numpy array (float)
            B_l(r,r') array at the specified value of l
        """
        s2 = padf.shape[2]
        z = np.cos( 2.0*np.pi*np.arange(s2)/float(s2) )
        Pn = legendre( int(l) )
        if self.legendre_norm:
            p = Pn(z)*np.sqrt(2*l+1)
        else:
            p = Pn(z)
        #print "pshape", p.shape, padf.shape

        mat = np.zeros( padf.shape )
        for i in np.arange(padf.shape[0]):
             for j in np.arange(padf.shape[1]):
                  mat[i,j,:] = p[:]

        output = np.sum( mat*padf, 2 )
        return output


    def project_padf_Legendre_zsamp( self, padf, l ):
        """        
        Compute B_l(r,r') matrix from padf for a particular value of l
        Uses regular cos(theta) sampling.

        Parameters
        ---------
        padf : vol object
            input padf volume

        l : int
            spherical harmonic order

        Returns
        -------
        output : numpy array (float)
            B_l(r,r') array at the specified value of l
        """
        s2 = padf.shape[2]
        
        z = np.arange(s2)*2/s2 - 1
        acs = (self.nth/2)*np.arccos(z)/np.pi
        acs2 = (self.nth/2)*np.arccos(z[::-1])/np.pi + self.nth/2
        zorder = 1 

        zcorr = np.zeros(padf.shape )
        for ic in np.arange(padf.shape[0]):
            for jc in np.arange(padf.shape[0]):
                    zcorr[ic,jc,:] = map_coordinates(padf[ic,jc,:], [acs], order=zorder) 
                    zcorr[ic,jc,:] +=  map_coordinates(padf[ic,jc,:], [acs2], order=zorder) 

        #z = np.cos( 2.0*np.pi*np.arange(s2)/float(s2) )
        Pn = legendre( int(l) )
        if self.legendre_norm:
            p = Pn(z)*np.sqrt(2*l+1)
        else:
            p = Pn(z)      
        #print "pshape", p.shape, padf.shape

        mat = np.zeros( padf.shape )
        for i in np.arange(padf.shape[0]):
             for j in np.arange(padf.shape[1]):
                  mat[i,j,:] = p[:]

        output = np.sum( mat*zcorr, 2 )
        return output




    def calcBlrr( self, padf, lmax ):
      """
        Computes all B_l(r,r') matrices from the padf volume
        
        Parameters
        ---------
        padf : vol object
            padf volume
    
        lmax : int
            maximum spherical harmonic order to compute

        Returns
        -------
        blrr : blqq object
            B_l(r,r') arrays for all values of l
      """
      blrr = []
      for l in np.arange(lmax):
          blrr.append( self.project_padf_Legendre(padf,l))
      return blrr


    #
    # apply the Ewald filter to a set of blqq matrices
    #
    """   def Blqq_Ewald_filter( self, blqqin, corrshape ):

      blqqout = blqq.blqqarray(self.nl, qmax=self.qmax, rmax=self.rmax,
                                     tablenzero=self.tablenzero,tablelmax=self.tablelmax) 
      
      lmax = blqqin.nl
      nq = blqqin.l[0].nq
      for iq in np.arange(nq):
         for iq2 in np.arange(nq):

             
              corr = np.zeros( corrshape[2] )
              qln1 = blqqin.jnuzeros[0,iq+1] / (2*np.pi*self.rmax)
              qln2 = blqqin.jnuzeros[0,iq2+1] / (2*np.pi*self.rmax)
                  
              for l in np.arange(2,lmax):

                   if (l%2)!=0:
                        continue

                   s2 = corrshape[2]
                   z = np.cos( 2.0*np.pi*np.arange(s2)/float(s2) )
                   Pn = legendre( int(l) )
                   p = Pn(z)
                               
                   corr[:] += blqqin.l[l].data[iq,iq2]*p[:]
    
              tmp = self.Blqq_legendre( qln1, qln2, corr, True )
              for l in np.arange(lmax//2):
                  blqqout.l[2*l].data[iq,iq2] = tmp[l]

      return blqqout
    """
    #
    # Generates self.blqq matrices from the correlation function
    # Each matrix is sampled at the zeros the corresponding spherical bessel function
    #
    def Blqq_calc_fast_interp(self, method='svd', order=3):
        """
        Compute the B_l(q,q') arrays from the correlation volume using radial q interpolation.

        The q-sampling for the l=0 spherical Bessel zeros are used.
        This 'fast' because the same q-sampling is used for each order l.        

        Parameters
        ---------
        method : str
            matrix inversion method for obtaining B_l matrices from correlation volume
            'svd' - singular value decomposition
            'legendre' - use orthogonality properties of Legendre polynomials
        
        order : int
            order to use for spline interpolation    

        Returns
        -------
            N/A
            Result stored in self.blqq and self.blqqtmp attributes
        

        """
        
        self.blqqtmp = blqq.blqqarray(self.nl, nq=self.blqq.l[0].nq, 
                                      tablenzero=self.tablenzero,tablelmax=self.tablelmax) 
        
        # define a list of matrices to store the Blqq matrices
        rsinvlist = []
        
        ###plt.figure()
        for l in np.arange(self.nl):
            mat = self.blqq.resampling_matrix(l, 0, self.blqq.l[l].nq, self.blqq.l[0].nq,\
                                                 self.qmax, self.rmax )
            
            u, s, vh = np.linalg.svd( mat, full_matrices=False )
            #print( "resamp smax smin", s[0], s[-1], "; dims u s vh", u.shape, s.shape, vh.shape )
            igood = np.where(s > 0.5)
            sinv = s*0.0
            sinv[igood] = 1.0/s[igood]
            inv = np.dot(np.dot(vh.transpose(), np.diag(sinv)), u.transpose())
            #print("resamp inv shape", inv.shape, np.max(inv), np.min(inv) )
#            inv = sp.linalg.pinv(mat)    # I MAY NOT HAVE TO INVER THIS, BUT JUST SOLVE STUFF
            rsinvlist.append(inv)

            ###plt.plot(s)
        ###plt.draw()
        ###plt.show()
            
        nq = self.blqq.l[0].nq
        #print("DEBUG BLQQ_CALC_FAST nq", nq)
        qlist = []
        qindices = np.zeros( (nq,nq,2) ).astype(np.int) - np.int(1)
        for iq in np.arange(nq):
           for iq2 in np.arange(nq):

              qln1 = self.blqq.jnuzeros[0,iq+1] / (2*np.pi*self.rmax)
              qln2 = self.blqq.jnuzeros[0,iq2+1] / (2*np.pi*self.rmax)
              #qln1 = self.blqq.jnuzeros[0,iq] / (2*np.pi*self.rmax)
              #qln2 = self.blqq.jnuzeros[0,iq2] / (2*np.pi*self.rmax)
              
              ic2 =  (self.nq)*(qln1/self.qmax) 
              jc2 =  (self.nq)*(qln2/self.qmax) 
              qlist.append( np.array([ic2, jc2]) )


              ic = np.round( self.nq*(qln1/self.qmax) ).astype(np.int)
              jc = np.round( self.nq*(qln2/self.qmax) ).astype(np.int)
              if (ic>=nq) or (jc>=nq): continue 
              qindices[iq,iq2,0] = ic
              qindices[iq,iq2,1] = jc             

        #print("DEBUG interpolating correlation onto legendre zeros")
        qlist = np.array(qlist) 
        #print( "qlist shape", qlist.shape )
        corr_interp = np.zeros( (nq, nq, self.nth) )
        corr_interp2 = np.zeros( (nq, nq, self.nth) )
        #print( self.corrvol.shape)
        for ith in np.arange( self.nth ):
             #print("\r ith", ith)
             corr_interp[:,:,ith] = map_coordinates( self.corrvol[:,:,ith], qlist.transpose(), order=order ).reshape( (nq,nq) )
         
        for iq in np.arange(nq):
           for iq2 in np.arange(nq):
                 if (qindices[iq,iq2,0]>=0) and (qindices[iq,iq2,1]>=0):
                      corr_interp2[iq,iq2,:] = self.corrvol[qindices[iq,iq2,0],qindices[iq,iq2,1],:]


        #print("qlist", qlist)

        #fname = "corr_interp_DEBUG.npy"
        #np.save(fname, corr_interp)
        #fname = "corr_interp2_DEBUG.npy"
        #np.save(fname, corr_interp2)

        #corr_interp = corr_interp2 #HACK


        #print("DEBUG BLQQ_CALC_FAST nq", nq)
        for iq in np.arange(nq):
           for iq2 in np.arange(nq):

              
              qln1 = self.blqq.jnuzeros[0,iq+1] / (2*np.pi*self.rmax)
              qln2 = self.blqq.jnuzeros[0,iq2+1] / (2*np.pi*self.rmax)
              
              # CALCULATE AND INVERT FMAT; MAY NOT HAVE TO INVERT JUST SOLVE
              if method=='svd':
                 #mat = self.fmat( qln1, qln2 )
                 #matinv = sp.linalg.pinv2( mat, cond=0.5)   
                 #tmp = np.dot( matinv.transpose(), corr_interp[iq,iq2,:])
             

                 zmat = self.fmat_zsamp( qln1, qln2 )
                 z = np.arange(self.nth)*2/self.nth - 1
                 acs = (self.nth/2)*np.arccos(z)/np.pi
                 acs2 = (self.nth/2)*np.arccos(z[::-1])/np.pi + self.nth/2
                 zorder = 1
                 zcorr = map_coordinates(corr_interp[iq,iq2,:], [acs], order=zorder) 
                 zcorr +=  map_coordinates(corr_interp[iq,iq2,:], [acs2], order=zorder) 

                 zmatinv = sp.linalg.pinv( zmat, cond=0.5 )
                 tmp = np.dot( zmatinv.transpose(), zcorr )
                 
 
              elif method=='legendre':
                 tmp = self.Blqq_legendre( qln1, qln2, corr_interp[iq,iq2,:], True )

              else:            
                 print("invalid method option for calculating Blqq. 'svd' or 'legendre'. Exiting.")
                 exit()
 
              for k in np.arange(self.nl//2):
                 self.blqqtmp.l[2*k].data[iq,iq2] = tmp[k]
                 
              for k in np.arange(self.nl//2-1):
                 self.blqqtmp.l[2*k+1].data[iq,iq2] = 0.0

        # RESAMPLE USING THE RESAMPLING MATRICES...
        for k in np.arange(self.nl//2):
           tmp = np.dot(rsinvlist[2*k],self.blqqtmp.l[2*k].data)
           self.blqq.l[2*k].data = np.dot( tmp, rsinvlist[2*k].transpose())

