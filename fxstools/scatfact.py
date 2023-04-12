"""
scatfact.py

Generate scattering factors for atomic elements

This script generates the values of atomic scattering factors on 1D and 2D
q-space grids

This file can not be executed independently. It is imported as a module 
and contains the following classes:

    * sf  - represents scattering factor of an atom
    * sfdata - scattering factor data for a group of elements in a sample
"""
import numpy as np
import os


class sf:
    """
    A class to represent the scattering factor of an atom
    
    Attributes
    ----------
    nx : int
        number of elements in 1D samples and side length of 2D samples
    sf1d : numpy array of size nx (floats)
        scattering factor as a function of abs(q) (=magnitude of the 
                                                   scattering vector q)
    sf2d : numpy array of size nx x nx (floats)
        scattering factor as a function of q on a 2D grid
    """
    def __init__(self,nx=100):
            """
            Constructor for the scattering factor of an atom
            
            Attributes
            ----------
            nx : int
                number of elements in 1D samples and side length of 2D samples
            sf1d : numpy array of size nx (floats)
                scattering factor as a function of abs(q) (=magnitude of the 
                                                           scattering vector q)
            sf2d : numpy array of size nx x nx (floats)
                scattering factor as a function of q on a 2D grid
            """
            self.nx = nx
            self.sf2d = np.zeros((nx,nx))
            self.sf1d = np.zeros(nx)



class sfdata:
    """
    scattering factor data for a group of elements in a sample
    
    Attributes
    ----------
    nx : int
        number of elements in 1D samples and side length of 2D samples
    
    wl : float
        wavelength of the radiation in metres (m)

    dz : float
        distance between the samplea and the detector in metres (m)

    pw : float
        width of a detector pixel in metres (m)
        
    cenx : float
        x-coordinate of the diffraction pattern centre in pixels. 
        Non-integer values are handled with interpolation of the pattern
   
    ceny : float
        y-coordinate of the diffraction pattern centre in pixels. 
        Non-integer values are handled with interpolation of the pattern
        
    henkeflag : Bool
        Sets whether Henke parameters are used for wavelength dependence
        of the scattering factor
        
    sfsource : str
        filename containing the scattering factor parameters
        
    sfparams : numpy array float
        array of scattering factor parameters read from file
        
    x : numpy array float
        x coordinates of detector pixels
        
    y : numpy array float
        y coordinates of detector pixels
    
    z : numpy array float
        z coordinates of detector pixels
    
    len : numpy array float
        realspace distance between sample origin and each detector pixel
        
    qx, qy, qz : numpy arrays float
        reciprocal space qx (qy, qz) coordinates for each detector pixel
    
    q2 : numpy array float
        square magnitude of each reciprocal vector (2D)
    
    qlen : numpy array float
        magnitude of each reciprocal vector (2D)
        
    q1D  : numpy array float
        magnitude of 1D regularly sampled q values
    
    q1Dsq : numpy array float
        q1D squared
        
    qindices : tuples
        lists of pixel indices (2D) for each reciprocal space magnitute (q)
    
    sflist : list 
        list of scattering factors for elements in the sample
        
    hdata : numpy array
        Henke data on wavelength dependence of the scattering factor
        
    Methods
    ---------
    
        
    
    """
    def __init__(self,nx=100,wl=1e-10,dz=0.1,pw=1e-5,
                 cenx=-1,ceny=-1,henkeflag=False,sfsource="Bwxray.fac"):
        """
        Constructor for scattering factor data
        
        Attributes
        ----------
        nx : int
            number of elements in 1D samples and side length of 2D samples
        
        wl : float
            wavelength of the radiation in metres (m)
        
        dz : float
            distance between the samplea and the detector in metres (m)
        
        pw : float
            width of a detector pixel in metres (m)
            
        cenx : float
            x-coordinate of the diffraction pattern centre in pixels. 
            Non-integer values are handled with interpolation of the pattern
        
        ceny : float
            y-coordinate of the diffraction pattern centre in pixels. 
            Non-integer values are handled with interpolation of the pattern
            
        henkeflag : Bool
            Sets whether Henke parameters are used for wavelength dependence
            of the scattering factor
            
        sfsource : str
            filename containing the scattering factor parameters
        """
    
        self.sfsource = sfsource
        self.nx = int(nx)
        self.wl = wl
        self.dz = dz
        self.pw = pw

        if cenx >= 0:
            self.cenx = cenx
        else:
            self.cenx = self.nx/2

        if ceny >= 0:
            self.ceny = ceny
        else:
            self.ceny = self.nx/2

        self.henkeflag = henkeflag

        #print("<scatfact.py>", self.nx, self.cenx, self.ceny )

        #
        # load the parameters
        #
        self.load_sf_parameters()
        self.calculate_qgrids()


    def load_sf_parameters(self):
        """Load the atomic scattering factor parameters from file
    
        Attributes required
        ----------
        self.sfsource - file name containing scattering factor parameters
    
        Attributes overwritten
        ----------
        self.sfparams = a numpy array containing the scattering factor parameters
        """
        dir = os.path.dirname(os.path.abspath(__file__))
        self.sfparams = np.loadtxt( dir+"/"+self.sfsource, delimiter=",")
        


    def calculate_qgrids( self ):
        """Calculate arrays indicating q-space coordinates of the detector pixels
    
        Attributes required
        ----------
        self.nx : int 
            number of pixels on size length of detector array (assumed square)

        self.pw : float
            pixel width (metres)

        self.cenx : float
            x coordinate of beam centre in pixels    

        self.ceny : float
            y coordinate of beam centre in pixels

        self.dz : float
            sample to detector distance (metres)

        self.wl : float
            x-ray wavelength (metres)   
 
        Attributes overwritten
        ----------
        self.x, self.y, self.z : numpy arrays, float
            pixel coordinates (metres)

        self.len : numpy array, float
            distance of pixel from beam centre (metres)

        self.qx, self.qy, self.qz : numpy arrays, float
            q-space coordinates of each detector pixel
            accouting for the Ewald sphere

        self.qlen, self.q2 : numpy arrays, float
            q-space distance of each pixel from beam centre (qlen)
            q-space distance squared (q2)


        self.q1d, self.q1dsq : numpy arrays, float
            q-space distance of each radial bin (q1d)
            q-space distance of each radial bin squared (q1dsq)

        self.qindices : list of tuples
            list of pixel indices for each qbin
        """
        xy = np.mgrid[:self.nx,:self.nx]
        #print(xy.shape, self.nx
        self.x = (xy[0] - self.cenx)*self.pw
        self.y = (xy[1] - self.ceny)*self.pw
        self.z = np.ones( (int(self.nx),int(self.nx)) )*self.dz
        self.len = np.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)
        lwl = self.len * self.wl*1e10

        self.qx = self.x / lwl
        self.qy = self.y / lwl
        self.qz =  - 1/(self.wl*1e10) + self.z/lwl

        self.q2 = self.qx*self.qx + self.qy*self.qy + self.qz*self.qz
        self.qlen = np.sqrt( self.q2 )
        
        self.q1d = np.max(self.qx) * np.arange(self.nx)/self.nx
        self.q1dsq = self.q1d**2
        
        self.qindices = []
        for iq in np.arange(self.nx-1):    
            self.qindices.append( np.where( (self.qlen>self.q1d[iq])*(self.qlen<self.q1d[iq+1]) ))

        
    # calculate the sf from a specific element
    def sf_calc( self, Z , elem ):
        """Calculate the atomic scattering factor for an element
           on the detector pixel array and the 1D radial bins 
    
        Parameters
        ----------
        Z : int
            atomic number of the element

        elem : string
            string that names Henke data file name without extension
            (extension assumed to be .nff) 
    
        Attributes required
        ----------
        self.sfparams : numpy array
             a numpy array containing the scattering factor parameters
            
        self.q2, self.q1dsq : numpy array, floats
            arrays of q-space distances squared. See self.calculated_qrids()

        self.nx : int
            number of pixels on size length of detector array

        self.wl : float
            wavelength of x-ray beam

        Returns
        ---------
        sfz : sf class
            contains scattering factors evaluated on 2D detector pixels
            and 1D radial bins.
        """
        sfz = sf( self.nx )
        if Z==0:
            if elem=="PP":
                p = np.zeros( 10 )
                p[0] = 1.0
                #print("debug point scatterer found.")
        else:
            iz = np.where( self.sfparams[:,0]==Z)[0][0]
            p = self.sfparams[iz,1:]
            #print("<scat_fact.py>", Z, elem, p, np.max(self.q2), np.min(self.q2) )
            

        if self.henkeflag:
            energy = 1.23984198e-6 / self.wl
            self.get_henke_data( elem )
            fvals = self.get_henke_f_vals( energy )

        for ic in np.arange(5):
            sfz.sf2d +=  p[ic]*np.exp( -p[ic+5]*self.q2/4)
    
            sfz.sf1d +=  p[ic]*np.exp( -p[ic+5]*self.q1dsq/4)

        if self.henkeflag:
            sfz.sf2d += fvals[0] - np.sum(p)
            sfz.sf1d += fvals[0] - np.sum(p)
            
        return sfz

    def sf_list_calc( self, zlist, elist ):
        """Generates scattering factors evaluated on the detector grid
           and 1D radial bins for all relevant atomic elements
    
        Parameters
        ----------
        zlist : list of ints
            atomic numbers of the elements in the sample

        elist : list of strings
            string that names Henke data file name without extension
            for all elements in the sample

        Attributes overwritten
        ----------
        self.sflist - list of sf class objects
            stores scattering factors for all elements in the sample
            evaluated on detector array and 1D bins
        """
        self.sflist = []

        for (Z, elem) in zip(zlist, elist):
            sfz = self.sf_calc(Z, elem)
            self.sflist.append( sfz )

                             
    def get_henke_data( self, elem ):
        """reads Henke data from file for a specific element
           required for wavelength dependent scattering factors
    
        Parameters
        ----------
        elem : string
            file name of Henke data file without extension (assumed .nff)
      

        Attributes overwritten
        ----------
        self.hdata - numpy array
            energies and Henke factors
        """
        fname = elem+".nff"
        hpath = os.path.dirname(os.path.abspath(__file__))+"/henke/"
        self.hdata = np.loadtxt( hpath+fname, skiprows=1 )

    # get a specific value for the wavelength being used
    def get_henke_f_vals( self, energy ): 
        """Computes wavelength dependent correction factors f0, f1
            using Henke data
        Parameters
        ----------
        energy : float
            photon energy in eV
 
        Returns
        ----------
        [f0, f1] : numpy array, floats
            f0 and f1 values
        """

        for i in np.arange(self.hdata.shape[0]-1):
            
            if (energy>self.hdata[i,0])and(energy<self.hdata[i+1,0]):
                factor = (energy - self.hdata[i,0])/ (self.hdata[i+1,0]-self.hdata[i,0]) 
                f0 = self.hdata[i,1] + factor*(self.hdata[i+1,1]-self.hdata[i,1]) 
                f1 = self.hdata[i,2] + factor*(self.hdata[i+1,2]-self.hdata[i,2]) 
                break

        return np.array([f0,f1])

    def sf1d_to_sf2d( self, sf1d ):
        """Maps the 1D scattering factor onto the 2D detector grid
           resulting in a radially symmetric 2D scattering factor.
    
        Parameters
        ----------
        sf1d : numpy array
            scattering factors evaluated on 1D radial bins

        Returns
        ---------
        sf2d : 2D numpy array
            scattering factors evaluated on 2D detector pixels

        Attributes required
        ----------
        self.nx : int
            side length of 2D detector array in pixels

        self.qindices : list of tuples
            indices of 2D pixels corresponding to each radial bin

        """
        sf2d = np.zeros( (self.nx,self.nx) )
        for iq, qind in enumerate(self.qindices):
            sf2d[qind] = sf1d[iq]
        return sf2d
