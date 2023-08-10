"""

correlationTools.py

Analysis scripts for performing angular correlation analysis
         
Classes
    angular_correlation
"""

import numpy as np
import scipy.ndimage as sdn


class angular_correlation:
    """Contains useful methods for calculating angular correlations
    """
    
    def polar_plot( self,data, nr, nth, rmin, rmax, thmin, thmax, cenx, ceny, submean=False ):
        """
	    Convert a 2D image into a r v theta plot

        Parameters
	----------
        data : numpy array
	    2D array of intensity data
	   
	nr : int
	    number of radial bins

	nth : int
            number of angular bins

	rmin : float
	    value of smallest radial bin. (arbitrary units)

	rmax : float
	    value of largest radial bin

	thmin : float
	    value of smallest angular bin (radians)

	thmax : float
	    value of largest angular bin (radians)

	cenx : float
	    row of beam centre (pixel units)

	ceny : float
	    column of beam centre (pixel units)

	submean : bool
	    subtract mean value in each radial ring if True

	Returns
	-------
	out : numpy array (float)
	    Data interpolated onto a r vs theta grid

	    
        """
        # r and theta arrays
        rarr = np.outer( np.arange(nr)*(rmax-rmin)/float(nr) + rmin, np.ones(nth) )
        tharr = np.outer( np.ones(nr), np.arange(nth)*(thmax-thmin)/float(nth) + thmin)
        
        newx = rarr*np.cos( tharr ) + cenx
        newy = rarr*np.sin( tharr ) + ceny
        
        newdata = sdn.map_coordinates( data, [newx.flatten(), newy.flatten()], order=3 )

        out = newdata.reshape( nr, nth )
        if submean == True:
            out = self.polar_plot_subtract_rmean( out  )

        return out


    def polar_plot_with_qbins( self,data, qbins, nth, thmin, thmax, cenx, ceny, submean=False, realout=True):
        """
        Convert a 2D image into a r v theta plot
        using irregular radial (q) bin locations

        Parameters
        ----------
        data : numpy array
            2D array of intensity data
      
        qbins : numpy array (float)
            coordinates of radial bins
 
        nth : int
            number of angular bins

        thmin : float
    '       value of smallest angular bin (radians)

        thmax : float
            value of largest angular bin (radians)

        cenx : float
            row of beam centre (pixel units)

        ceny : float
            column of beam centre (pixel units)

        submean : bool
            subtract mean value in each radial ring if True

        Returns
        -------
        out : numpy array (float)
            Data interpolated onto a r vs theta grid
        """

        nr = qbins.size
        #print( "pp debug", data.shape, qbins.size, nth ) 

        # r and theta arrays
        rarr = np.outer( qbins, np.ones(nth) )
        tharr = np.outer( np.ones(nr), np.arange(nth)*(thmax-thmin)/float(nth) + thmin)
        
        newx = rarr*np.cos( tharr ) + cenx
        newy = rarr*np.sin( tharr ) + ceny
        #print( "pp debug newx newy", np.min(newx), np.max(newx), newx.shape)
        np.save( "newx.npy", newx)
        np.save( "newy.npy", newy)
        newdata = sdn.map_coordinates( data, [newx.flatten(), newy.flatten()], order=3 )

        out = newdata.reshape( nr, nth )
        if submean == True:
            out = self.polar_plot_subtract_rmean( out  )

        if realout: out = np.real(out) 
        return out


    def qbins( self, nq, pmax, dz, wl, pw ):
        """
        Generates a list of q bins based on Ewald sphere curvature

        Parameters
        ----------
        nq : int
            number of q bins

        pmax : int
            number of detector pixels from centre to detector edge

        dz : float
            sample to detector distance

        wl : float
            wavelength

        pw : float
            width of a detector pixel

        Returns
        -------
        qpix : numpy array (float)
            qvalues of each radial (q) bin

        
        """

        qmax = (2/wl)*np.sin( np.arctan(pmax*pw/dz)/2.0 )

        qind = np.arange(nq)*qmax/float(nq)

        qpix = (dz/pw)*np.tan(2.0 * np.arcsin( qind*(wl/2.0) ))

        return np.floor(qpix)


    def polar_plot_subtract_rmean( self, pplot ):
        """
        Subtract the mean value in each q-ring from a polar plot

        Parameters
        ----------
        pplot : numpy array (float)
            input polar plot

        Returns
        -------
        out : numpy array (float)
            polar plot with q-ring mean value subtracted
        """

        av = np.average( pplot, 1 )
        out = pplot -np.outer( av, np.ones( pplot.shape[1] ) )
        return out

    #
    # performs the angular correlation of each q-shell with itself
    #
    def polarplot_angular_correlation( self, polar, polar2=None):
        """
        Calculate the 2D angular correlation of a polar plot
        or cross-correlation of two polar plots

        Parameters
        ----------
        polar : numpy array (float)
            input polar plot

        polar2 : numpy array (float)
            second input polar plot. 
            If provided, then cross-correlation is computed
            between polar and polar2.
            If polar2 not provided, then auto-correlation of
            polar is computed

        Returns
        -------
        out : numpy array (float)
            angular correlation function        
        """

        fpolar = np.fft.fft( polar, axis=1 )

        if polar2 != None:
            fpolar2 = np.fft.fft( polar2, axis=1)
            out = np.fft.ifft( fpolar2.conjugate() * fpolar, axis=1 )
        else:
            out = np.fft.ifft( fpolar.conjugate() * fpolar, axis=1 )
       
        return out

    #
    # angular correlation of each q-shell with all other q-shells
    #
    def polarplot_angular_intershell_correlation( self, polar, polar2=None, realout=True):
        """
        Calculate the 3D angular correlation [C(q,q',theta)] of a polar plot
        or cross-correlation of two polar plots.
        Each q-ring is correlated with every other q-ring.

        Parameters
        ----------
        polar : numpy array (float)
            input polar plot

        polar2 : numpy array (float)
            second input polar plot. 
            If provided, then cross-correlation is computed
            between polar and polar2.
            If polar2 not provided, then auto-correlation of
            polar is computed

        realout : bool
            enure output if real valued if True.

        Returns
        -------
        out : numpy array (float)
            angular correlation function        
        """

        fpolar = np.fft.fft( polar, axis=1 )

        if np.any(polar2) != None:
            fpolar2 = np.fft.fft( polar2, axis=1)
        else:
            fpolar2 = fpolar
      
        out = np.zeros( (polar.shape[0],polar.shape[0],polar.shape[1]) , dtype=np.complex128)
        for i in np.arange(polar.shape[0]):
            for j in np.arange(polar.shape[0]):
                out[i,j,:] = fpolar[i,:]*fpolar2[j,:].conjugate()
        out = np.fft.ifft( out, axis=2 )

        if realout: out = np.real(out) 
        return out

        
    def apply_mask( self, func, mask ):
        """
        Multiplies an array by a mask array

        Parameters
        ----------
        func : numpy array
            numpy array of data

        mask : numpy array
            mask array containing 0s and 1s

        Returns
        -------
        func*mask : numpy array
        """
        return func*mask


    def mask_correction( self, corr, maskcorr ):
        """
        Corrected correlation function for effects of the mask.
        Divides corr by maskcorr wherever maskcorr is greater than 0.

        Parameters
        ----------
        corr : numpy array
            correlation function

        maskcorr : numpy array
            correlation of mask function

        Returns
        -------
        corr : numpy array
            correlation data divided by mask correlation
        """
        imask = np.where( maskcorr != 0 )
        corr[imask] *= 1.0/maskcorr[imask]
        return corr

    #    
    # pairwise correlation of (flattened) arrays
    #
    # not for angular correlations; good for correlation of mean asic values
    #
    def allpixel_correlation( self, arr1, arr2 ):
        """
        Returns the outer product between the flattened
        arr1 and arr2. 
        """
        out = np.outer( arr1.flatten(), arr2.flatten() )
        return out

    # pearson correlation of a 2D area
    def pearsonCorrelation2D( self, arr1, arr2, lim=None):
        """
        Computes the Pearson correlation between two numpy arrays

        Parameters
        ----------
        arr1, arr2 : numpy arrays (floats)
            Arrays with the same number of elenments

        Returns
        -------
        pc : float
            Pearson correlation value
        """
        if lim == None:
            lim = [0, arr1.shape[0], 0, arr1.shape[1]]
      
        a1 = arr1[lim[0]:lim[1],lim[2]:lim[3]]
        a2 = arr2[lim[0]:lim[1],lim[2]:lim[3]]
        
        c1 = a1 - np.average(a1)
        c2 = a2 - np.average(a2)
        pc = np.sum( c1*c2 ) /np.sqrt( np.sum(c1*c1) * np.sum(c2*c2))
        return pc

# returns pearson correlation of each q ring
    def pearsonCorrelation2D_angular( self, arr1, arr2, lim=None):
        """
        Computes the Pearson correlation between two polar plots
        as a function of radial q-bin. 

        Parameters
        ----------
        arr1, arr2 : numpy arrays (floats)
            Arrays with the same number of elenments

        Returns
        -------
        pc : numpy array
            Pearson correlation values as a function of q (1D)
        """
        if lim == None:
            lim = [0, arr1.shape[0], 0, arr1.shape[1]]
      
        a1 = arr1[lim[0]:lim[1],lim[2]:lim[3]]
        a2 = arr2[lim[0]:lim[1],lim[2]:lim[3]]
        
        c1 = a1 - np.outer( np.average(a1, 1), np.ones( a1.shape[1]) )
        c2 = a2 - np.outer( np.average(a2, 1), np.ones( a2.shape[1]) )
        pc = np.sum( c1*c2, 1 ) /np.sqrt( np.sum(c1*c1, 1) * np.sum(c2*c2, 1))
        return pc

        
