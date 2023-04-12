"""
thmask.py

Create angular masks to exclude certain angular ranges

Functions
    costheta_mask_sphB
    mask_costheta_mask
    mask_costheta_mask_v2

"""
import numpy as np
from scipy.special import legendre

def costheta_mask_sphB( lmax, costhmask ):
        """Make a mask for costheta sampling
            based on a finite number of spherical harmonics
            
        Parameters
        ----------
        lmax : int
            maximum spherical harmonic order l to include

        costhmask : numpy array
            mask function for theta
    
        Result
        ------
        Amatrix : numpy array (float)
            matrix to apply mask to spherical harmonic coefficients

        """

        nth = costhmask.size
        #z = np.cos( 2.0*np.pi*np.arange(nth)/float(nth) )
        z = 2*np.arange(nth)/float(nth) - 1
        Amatrix = np.zeros( (lmax,lmax) )
        for l in np.arange(lmax):
            for l2 in np.arange(lmax):

                Pn = legendre( int(l) )
                Pn2 = legendre( int(l2) )
                
                p = Pn(z)
                p2 = Pn2(z)

                Amatrix[l,l2] = np.sum( p*p2*costhmask ) * np.sqrt((2*l+1) * (2*l2+1)) / (float(nth))
        
        return Amatrix


def make_costheta_mask( nth, thmin, thmax ):
    """Make a cos(theta) mask based on theta bounds
       Note: Theta bounds must be between 0-180 degrees
            
    Parameters
    ----------
    nth : int
        number of theta samples

    thmin : float
        minimum theta value

    thmax : float
        maximum theta value

    Result
    ------
    costhmask : numpy array (float)
        mask that applies to theta-sampled arrays
    """
    fact = np.pi/180.0
    z = 2*np.arange(nth)/float(nth) - 1
    costhmask = np.ones(nth)    
    zmin =  np.cos(thmin*fact)
    zmax =  np.cos(thmax*fact)

    if zmin < zmax:
        costhmask[(z>zmin)*(z<zmax)] = 0.0     
    else: 
        costhmask[(z<zmin)*(z>zmax)] = 0.0     
    
    return costhmask

def make_costheta_mask_v2( nth, zmin, zmax ):
    """Make a z=cos(theta) mask based on z  bounds
       Note: Theta bounds must be between - 1 and 1
            
    Parameters
    ----------
    zmin : int
        minimum z=cos(theta) value

    zmax : float
        maximum z=cos(theta) value

    Result
    ------
    costhmask : numpy array (float)
        mask that applies to z-sampled arrays
    """
    z = 2*np.arange(nth)/float(nth) - 1
    costhmask = np.ones(nth)    

    if zmin < zmax:
        costhmask[(z>zmin)*(z<zmax)] = 0.0     
    else: 
        costhmask[(z<zmin)*(z>zmax)] = 0.0     
    return costhmask


if __name__ == '__main__':

    nth = 1000
    lmax = 64
    #costhmask = np.ones(nth)
    costhmask = make_costheta_mask( nth, 0, 20 )
    Amatrix = costheta_mask_sphB( lmax, costhmask )
    for l in np.arange(lmax):
        print( l, Amatrix[l,l] )    
