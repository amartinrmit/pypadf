"""
vol.py

Store data in n dimensions (e.g. spatial) with dimension names and range limits.
"""


import numpy as np


class dimension:
    """
    Stores information about a 1D spatial dimension
    
    Attributes
    ----------
    name : string
        name of the dimension, e.g 'x'
        
    npix : int
        number of sampling point in the dimension

    dmin : float
        coordinate value of the first sampling point

    dmax : float
        coordiante value of the last sampling point

    """

    def __init__(self,npix2=10,dmin2=0,dmax2=9,name="dim"):
        """Constructs dimension class and sets parameter values
        """

        self._name = name
        self._npix = npix2
        self._dmin = dmin2
        self._dmax = dmax2

    def pnts(self):
        """generates array of coordinate values for a dimension
        
         Returns:
            p : numpy array of floats 
                coordinate values evenly spaced from dmin to dmax         
        """
        p = (np.arange(self.n)*(self.max - self.min)/self.n) + self.min
        return p


    """
    @property
    def name(self):
        return self._name


    @property
    def npix(self):
        return self._npix

    @property
    def dmin(self):
        return self._dmin

    @property
    def dmax(self):
        return self._dmax
    """        
    
    

class Vol2:
    """Stores a data volume in n dimensions

    Attributes:
    ----------
    dimnames : list of strings
        name of each dimension

    dimlen : list of ints
        number of sampling points for each dimension

    dmin : list of floats
        coordinate value of the first sampling point for each dimension

    dmax : list of floats
        coordiante value of the last sampling point for each dimension

    dims : list of dimension objects
        stores information about each dimension

    vol : len(dimnames) dimensional numpy array
        stores data volume
    """

    def __init__(self,dimnames=None,dimlen=None,dmin=None,dmax=None):
        """ Construct the Vol2 class
                Dimension class created for each dimension
                Volume attribute (n dimension numpy array) created
        """
        self.ndims = len(dimnames)
        if (len(dimlen)!=self.ndims) or (len(dimnames)!=self.ndims) \
                or (len(dmin)!=self.ndims) or (len(dmax)!=self.ndims):
            print( "incorrect initialization of Vol - all arguments must be of same length")
            exit()

        self.dims = {}
        for i in np.arange(self.ndims):
            d = dimension( dimlen[i], dmin[i], dmax[i], dimnames[i] )
            self.dims.update( {d._name : d} ) 

        self.vol = np.zeros( dimlen ) 
        
