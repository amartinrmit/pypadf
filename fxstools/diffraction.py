"""
diffraction.py

Tools to calculate a diffraction pattern from a pdb file

Classes:
    diffraction

Functions:
    fast_diffraction
    atom_diffraction
    vec_norm
    
Classes:
    unit_cell
    diffraction
"""
import numpy as np
import matplotlib.pyplot as plt
import fxstools.pypdb as pypdb
import fxstools.scatfact as sf
from fxstools.scatfact import sfdata
from fxstools.quaternions import rotation_matrix
from numba import jit
import time
#import jax

@jit(nopython=True, cache=True)
def fast_diffraction( nx, alist, qx, qy, qz, tmpr, tmpi ):

    for v in alist:
        
        # add option to rotate the vector
        #        if self.rotflag:
        #            v = np.dot( rmat, v )
            
        # dot product
        dotp = v[0]*qx + v[1]*qy + v[2]*qz
        tmpr = tmpr + np.cos( 2*np.pi*dotp )
        tmpi = tmpi + np.sin( 2*np.pi*dotp )
    return [tmpr, tmpi]


@jit(nopython=True, cache=True)
def atom_diffraction( v, qx, qy, qz ):
    """compute expontial for far-field diffracction with numba
    
    Parameters
    ----------
    v : numpy array of floats
        an atom's 3D position vector. 1D array of 3 floats.
    
    qx, qy, qz : float
        q-space coordinates of a reciprocal lattice vector
    """
    # dot product
    dotp = v[0]*qx + v[1]*qy + v[2]*qz
    #        tmpr = tmpr + np.cos( 2*np.pi*dotp )
    #        tmpi = tmpi + np.sin( 2*np.pi*dotp )
    return np.exp( 2 * np.pi * 1j * dotp)

@jit(nopython=True, cache=True)
def vec_norm( v ):
    """length of a vector (using numba)
    """
    return np.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])


class unitcell():
    """
    Contains unit cell parameters

    Parameters
    ----------
    
    a, b, c : float
        lengths of lattice vectors a, b, c

    alpha, beta, gamma : float
        angles between lattice vectors (following standard naming conventions)
    """
    def __init__(self,a=-1.0,b=-1.0,c=-1.0,alpha=90,beta=90,gamma=90):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma    



class diffraction(sfdata):
    """
    Tools to compute a diffraction pattern from a pdb file 

    Parameters
    ----------

    pdbname : str
        name of the .pdb file

    outpath : str
        location to save output files

    tag : str
        text to prefix output files (change this to avoid overwriting old simulations)

    fext : str
        file extension for output files (e.g. .npy)

    nx : int
        side length of detector in pixels (assumed square array)

    wl : float
        wavelength of incident beam in metres

    dz : float
        sample to detector distance (metres)

    pw : float
        physical width of a detector pixel (metres)

    cenx : int
        x (row) coordinate of the beam centre

    ceny : int
        y (column) coordinate of the beam centre

    henkeflag : Bool
        Use Henke parameters for wavelength dependent scattering factors

    axis : numpy array (floats)
        Axis of sample rotation. 3 element numpy array.

    theta : float
        rotation angle in radians

    rotflag : Bool
        Set to rotate the sample

    rmax : float
        Exclude all atoms outside a radius of rmax from the origin
        (Units are same as coordinate units in the pdb file)

    npulse : float
        estimated number of photons per exposure (per pulse)
        
    beamarea : float
        estimated beam area in m^2
    """
    def __init__(self,pdbname,outpath=None,tag=None,fext=None,
                 nx=100,wl=1e-10,dz=0.1,pw=1e-5,
                 cenx=-1,ceny=-1,henkeflag=False,axis=np.array([1.0,0.0,0.0]),
                 theta=0.0,rotflag=False,rmax=-1, npulse=1e11, beamarea=1e-14):
        
        self.pdbname = pdbname
        self.outpath = outpath
        self.tag     = tag
        self.fext    = fext
        sfdata.__init__(self,nx=nx,wl=wl,dz=dz,pw=pw,cenx=cenx,ceny=ceny,henkeflag=henkeflag)

        #
        # load the pdb file  
        #   
        self.load_pdb()

        self.rotflag = rotflag
        self.axis = axis
        self.theta = theta
        self.rmax = rmax
        self.beamarea = beamarea
        self.npulse = npulse

        self.classical_erad = 2.8179403267e-15

        np.random.seed(None)

    def load_pdb(self):
        """Read the atom positions from the .pdb file
        """
        #
        # load the structure from the pdb file
        #
        if self.pdbname[-4:]=='.pdb':
            self.pdb = pypdb.pdb(self.pdbname)
            self.pdb.read_pdb()
        else:
            self.blah = np.loadtxt(self.pdb)
        #
        # calculate element scattering factors
        #
        sfdata.sf_list_calc(self, self.pdb.zlist, self.pdb.elements )

    def diffraction2D(self):
        """Calculate the 2D diffraction pattern
        """
        #intialise scattered wave
        self.swave2d = np.zeros( (self.nx,self.nx), dtype=np.complex128 )

        if self.rotflag:
            rmat = rotation_matrix( self.axis, self.theta )

        for ie in np.arange(len(self.pdb.elements)):

 
            """
            tmp = np.zeros( (self.nx,self.nx), dtype=np.complex128 )
            for a in self.pdb.sorted_atom_list[ie]:
                

                
                v = np.array([a.x, a.y, a.z])

                # add option to rotate the vector
                if self.rotflag:
                    v = np.dot( rmat, v )

                # dot product
                dotp = v[0]*self.qx + v[1]*self.qy + v[2]*self.qz
                
                tmp += np.exp( 2*np.pi*1j*dotp )
            """

 
            """
            #
            # first jit attempt
            #
            alist = []
            for a in self.pdb.sorted_atom_list[ie]:
                alist.append( np.array([a.x, a.y, a.z]) )
        
            tmpr = np.zeros( (self.nx,self.nx) )
            tmpi = np.zeros( (self.nx,self.nx) )
            tmpr,tmpi = fast_diffraction(self.nx, alist,
                                   self.qx,self.qy,self.qz,tmpr,tmpi)
            tmp = tmpr + 1j*tmpi
            """        
            
            #
            # 2nd jit attempt
            #
            tmp = np.zeros( (self.nx,self.nx), dtype=np.complex128 )
            for a in self.pdb.sorted_atom_list[ie]:
                
                v = np.array([a.x, a.y, a.z])
                if (self.rmax>0)and(vec_norm(v)>self.rmax): continue

                if self.rotflag:
                    v = np.dot( rmat, v )

                # dot product
                               
                tmp += atom_diffraction( v, self.qx, self.qy, self.qz )*np.exp( - 0.25*a.Bfactor*self.q2 )
             
            #print("test", np.max(tmp), np.max(self.sflist[ie].sf2d), len(self.pdb.sorted_atom_list[ie]))    
            self.swave2d += tmp*self.sflist[ie].sf2d
            

        self.dp2d = np.abs(self.swave2d)**2
        self.dp2d *= self.npulse*(self.classical_erad**2)/self.beamarea
        self.dp2d = self.solid_angle_correction(self.dp2d)
        #print("length pdb.elements", len(self.pdb.elements), np.max(self.dp2d))
        #print("<diffraction.py> - ", self.pdb.elements[ie], np.max(self.dp2d))



    def diffraction1D(self):
        """Calculate the 1D radial diffraction by integration of the 2D pattern
           (requires the 2D diffraction pattern to be calculated already)
        """        
        self.dp1d = np.zeros( self.nx )

        for i, iq in enumerate(self.qindices):
            if len(iq[0])>0:
                self.dp1d[i] = np.sum(self.dp2d[iq])/len(iq[0])

            
    def waxs2Dcalc(self):
        """Compute an isotropic 2D diffraction pattern from the 1D radial intensity
        """
        self.waxs2d = np.zeros( (self.nx,self.nx) )
        
        for iq in enumerate(self.qindices):
            self.waxs2d[iq] = self.dp1d[i] 
    
    def circ_shift_coordinates(self, uc, mins, recentre=False, rmax=-1):
        """Translate the atom positions with periodic boundar conditions

           Parameters
           ----------
           uc : unit_cell class instance
                unit cell parameters (can be an artificial box)

           mins : numpy array (floats)
                values of the minimum x,y,z coordinates in the unit cell

           recentre : Bool
                recentre the shifted unit cell

           rmax : float
                exclude atoms outside a radius of rmax from the origin
        """
        tmplist = []
        for i, atom in enumerate(self.pdb.atomlist):
            r = np.random.rand(3)
            if uc.a>0: self.pdb.atomlist[i].x = (self.pdb.atomlist[i].x + mins[0] + r[0]*uc.a)%uc.a - mins[0] 
            if uc.b>0: self.pdb.atomlist[i].y = (self.pdb.atomlist[i].y + mins[1] + r[1]*uc.b)%uc.b - mins[1]
            ## todo: if uc.c is 0 then we get a runtime error
            if uc.c>0: self.pdb.atomlist[i].z = (self.pdb.atomlist[i].z + mins[2] + r[2]*uc.c)%uc.c - mins[2]
            if recentre==True:
                self.pdb.atomlist[i].x += -uc.a/2
                self.pdb.atomlist[i].y += -uc.b/2
                self.pdb.atomlist[i].z += -uc.c/2
            if rmax>0:
               if np.sqrt(self.pdb.atomlist[i].x**2 + self.pdb.atomlist[i].y**2 + self.pdb.atomlist[i].z**2 )<rmax:
                    tmplist.append(self.pdb.atomlist[i])
            else: 
                    tmplist.append(self.pdb.atomlist[i])
        self.pdb.atomlist = tmplist 
        self.pdb.sort_atom_list()

    def solid_angle_correction(self,dp):
        """applies a correction for the solid angle of each pixel
           based on pixels oriented along the scattering vector"""
        xy = np.mgrid[0:self.nx, 0:self.nx]
        x = (xy[0]-self.cenx)*self.pw
        y = (xy[1]-self.ceny)*self.pw
        costh = np.cos( np.arctan( np.sqrt(x**2 + y**2)/self.dz))
        dpout = dp*costh*(self.pw/self.dz)**2 
        return dpout

    def polarisation_factor(self,vec):
        """applies a polarisation correction
           vec : numpy array; 2 floats
                direction of polarisation vector perpendicular to the beam
                first element is the horizontal direction; second is vertical
           dp : numpy array
                diffraction pattern to be corrected 
        """
        xy = np.mgrid[0:self.nx, 0:self.nx]
        x = (xy[0]-self.cenx)*self.pw
        y = (xy[1]-self.ceny)*self.pw
        theta = np.arctan( np.sqrt(x**2 + y**2)/self.dz)
        phi = np.angle( x+1j*y)
        pol = v[0]*(1-(np.sin(phi)*np.sin(theta))**2) + v[1]*(1-(np.cos(phi)*np.sin(theta))**2)
        return pol

    def poisson_sample(self, dp):
        return np.random.poisson(dp)

    def gaussian(self, nr, wid):
        r = np.arange(nr) - nr/2
        g = np.exp(- r*r/(2*wid*wid))
        g = np.roll(g, nr//2)
        g *= 1.0/np.sum(g)
        return g
  
    def convolve_gaussian( self, data, wid):
        g = self.gaussian( data.size, wid)
        fg = np.fft.fft(g)
        fdata = np.fft.fft(data)
        fout = np.fft.ifft( fdata*fg.conjugate())
        return fout



    def box_pair_distribution(self, nxin, rmax, nr, xlim, ylim, zlim):
        """Estimate the background contribution from the mean density to the radial distribution function
        """
        #print("box pair debug limits:", xlim, ylim, zlim) 
        nx, ny, nz = nxin, nxin, nxin
        xyz = np.mgrid[:nx,:ny,:nz]

        x = xyz[0]*xlim/nx #-xlim/2 
        y = xyz[1]*ylim/ny #-ylim/2
        z = xyz[2]*zlim/nz #-zlim/2
        rvol = np.sqrt( x*x + y*y + z*z)
        vol = np.zeros( (nx,ny,nz) )
        vol = (xlim - np.abs(x)) * (ylim - np.abs(y)) * (zlim - np.abs(z))

        # make the 1D r array and find the reverse indices.
        r = np.arange(nr)*rmax/nr
         
        pd, be = np.histogram( rvol.flatten(), bins=nr, range=(0,rmax), weights=vol.flatten() )
                
        """
        breakflag = False
        mind = np.min([xlim,ylim,zlim])
        limits = [2.0, 5.0, 10.0, 20.0,mind)
        for i, lim in enuerate(limits):
            if lim>mind:
                lim = mind
                breakflag = True
            if i==(len(limits)-1):
                tmp_xlim, tmp_ylim, tmp_zlim = xlim, ylim, zlim
            else:
                tmp_xlim, tmp_ylim, tmp_zlim = lim, lim, lim

            x = xyz[0]*tmp_xlim/nx #-xlim/2 
            y = xyz[1]*tmp_ylim/ny #-ylim/2
            z = xyz[2]*tmp_zlim/nz #-zlim/2
            rvol = np.sqrt( x*x + y*y + z*z)
            vol = np.zeros( (nx,ny,nz) )
            vol = (xlim - np.abs(x)) * (ylim - np.abs(y)) * (zlim - np.abs(z))

            # make the 1D r array and find the reverse indices.
            r = np.arange(nr)*rmax/nr
         
            pd, be = np.histogram( rvol.flatten(), bins=nr, range=(0,rmax), weights=vol.flatten() )
            if i==0:
                pdsum = pd
            else:
                pdsum += pd   ###FIX WEIGHTS BEFORE I USE THIS!!!!!

    
            if breakflag: break
        """
        return pd
             
        
    
    
    def saxs(self, nr=-1, box_nx=350, norm_rmin=10, norm_rmax=20):
        """Calculate the SAXS pattern from the pair distribution function
        """
        #intialise scattered wave
        nelem = len(self.pdb.elements)
        rmax = self.nx/self.q1d[-1]
        if nr<0: nr = self.nx*10
        r = np.arange(nr)*rmax/nr
        rinv = r*0.0
        rinv[1:] = 1.0/r[1:]

        self.partial_pd = np.zeros( (nelem,nelem,nr), dtype=np.complex128 )

        # put atoms into arrays
        xyz_arrays = []
        for ie in np.arange(len(self.pdb.elements)):
            xyz = np.zeros( (len(self.pdb.sorted_atom_list[ie]),3))
            for i, a in enumerate(self.pdb.sorted_atom_list[ie]):
                xyz[i,:] = np.array([a.x,a.y,a.z])
            xyz_arrays.append(xyz) 

        # compute the box pair distribution
        xmin, xmax, ymin, ymax, zmin, zmax = 1e6, -1e6, 1e6, -1e6, 1e6, -1e6
        for xyz in xyz_arrays:
            if np.max(xyz[:,0])>xmax: xmax=np.max(xyz[:,0])
            if np.max(xyz[:,0])<xmin: xmin=np.min(xyz[:,0])
            if np.max(xyz[:,1])>ymax: ymax=np.max(xyz[:,1])
            if np.max(xyz[:,1])<ymin: ymin=np.min(xyz[:,1])
            if np.max(xyz[:,2])>zmax: zmax=np.max(xyz[:,2])
            if np.max(xyz[:,2])<zmin: zmin=np.min(xyz[:,2])
        print("xyz min max", xmin, xmax, ymin, ymax, zmin, zmax)            
        xlim, ylim, zlim  = xmax-xmin, ymax-ymin, zmax-zmin

        print("computing box pair distribution (background)")
        t0 = time.perf_counter()
        self.pdbox = self.box_pair_distribution( box_nx, rmax, nr, xlim, ylim, zlim )
        self.pdbox[r>0] *= 1.0/(r[r>0]**2)
        print("done - time taken for box pair distribution:", time.perf_counter()-t0, "s" )

        #plt.figure("r")
        #plt.plot(r)
        #plt.draw()
        #plt.show()
        irnorm = (r>norm_rmin)*(r<norm_rmax) 

        pdboxnorm = np.sum(self.pdbox[irnorm])

        #
        # Compute the partial pair distributions
        #
        print("Computing the partial pair distributions...")
        start = time.time()
        for ie in np.arange(len(self.pdb.elements)):
            for ie2 in np.arange(len(self.pdb.elements)):
                if ie==ie2:
                    m=1
                elif ie2>ie:
                    m=2
                else:
                    continue
                    
                icount = 0  
                nchunkmax = 1000
                nchunk1 = np.min([nchunkmax,len(self.pdb.sorted_atom_list[ie])])
                nchunk2 = np.min([nchunkmax,len(self.pdb.sorted_atom_list[ie2])])
                niter = len(self.pdb.sorted_atom_list[ie])//nchunk1
                niter2 = len(self.pdb.sorted_atom_list[ie2])//nchunk2
                if niter%nchunkmax  != 0: niter+=1    
                if niter2%nchunkmax != 0: niter2+=1    

                for i in range(niter):
                    chunk1 = xyz_arrays[ie][i*nchunk1:(i+1)*nchunk1,:]
                    nc1 = chunk1.shape[0]
                    for j in range(niter2):
                        if (ie==2) and (ie2==2): print("ij", i, j)
                        if i==j:
                            m2 = 1
                        elif j<i:
                            m2 = 2
                        else:
                            continue

                        chunk2 = xyz_arrays[ie2][j*nchunk2:(j+1)*nchunk2,:]
                        nc2 = chunk2.shape[0]
                        x = (np.outer(chunk1[:,0], np.ones(nc2)) - np.outer(np.ones(nc1),chunk2[:,0])).flatten()
                        y = (np.outer(chunk1[:,1], np.ones(nc2)) - np.outer(np.ones(nc1),chunk2[:,1])).flatten()
                        z = (np.outer(chunk1[:,2], np.ones(nc2)) - np.outer(np.ones(nc1),chunk2[:,2])).flatten()
                        rchunk = np.sqrt(x*x + y*y + z*z)
                        h, be =  np.histogram( rchunk, bins=nr, range=(0,rmax)) 
                        if m==1:
                            self.partial_pd[ie,ie2,:] += m2*h
                        elif m==2:
                            self.partial_pd[ie,ie2,:] += m2*h
                            self.partial_pd[ie2,ie,:] += m2*h
        
        end = time.time()
        print("time to calculate pair distributions: ", end-start)

        #plt.figure()
        #plt.plot( r, self.partial_pd[0,0,:])
        #plt.title("pd")
        #plt.draw()
        #plt.show()
        print( "Peforming SAXS calculation")
        print( "elements used in saxs calculation", self.pdb.elements) 
        #
        # Compute the SAXS pattern
        #   
        self.saxs = np.zeros( self.nx ).astype(np.float64)
        nelem = len(self.pdb.elements)
        first_term = 0
        self.saxs_partial = np.zeros( (nelem,nelem,self.nx) ).astype(np.float64)
        self.normalisation = np.zeros( self.nx ).astype(np.float64)
        for ie in np.arange(len(self.pdb.elements)):
            sf = self.sflist[ie].sf1d**2
            self.saxs += np.real(len(self.pdb.sorted_atom_list[ie])*sf)
            first_term +=  np.real(len(self.pdb.sorted_atom_list[ie])*sf)
            self.normalisation += np.real(len(self.pdb.sorted_atom_list[ie])*sf)*(nr/rmax)
            for ie2 in np.arange(len(self.pdb.elements)):
                self.partial_pd[ie,ie2,r>0] *= 1.0/(r[r>0]**2)
                norm = np.sum(self.partial_pd[ie,ie2,irnorm])/pdboxnorm         
                self.partial_pd[ie,ie2,:] += -norm*self.pdbox
                self.partial_pd[ie,ie2,:] *= nr/rmax
                self.partial_pd[ie,ie2,:10] *= 0.0
                if ie==ie2:
                    m=1
                elif ie2<ie:
                    m=2
                else:
                    continue
                    
                sf = self.sflist[ie].sf1d*self.sflist[ie2].sf1d
                tmp = np.zeros( self.nx )
                for ir in np.arange(nr):
                   tmp += m*np.real(sf*self.partial_pd[ie,ie2,ir]*np.sin( 2*np.pi*self.q1d*r[ir])*r[ir])
                tmp[self.q1d>0] *= 1.0/self.q1d[self.q1d>0]
                self.saxs += tmp
                self.saxs_partial[ie,ie2,:] = tmp
        self.saxs *= 1.0/self.normalisation

        for ie in np.arange(len(self.pdb.elements)):
            for ie2 in np.arange(len(self.pdb.elements)):
                self.saxs_partial[ie,ie2,:] += first_term
                self.saxs_partial[ie,ie2,:] *= 1./self.normalisation


