"""
padfplot.py

Plotting tools for padf and correlation volumes

Classes:
    padfplot_dims
    
Functions:
    mult_radial_polynomial
    remove_angular_average
    multiply_by_sintheta
    generate_unit_labels
    make_gaussian_n
    convolve_gaussian_n
    corr_rescale
    extract_section
"""


import numpy as np
from scipy.ndimage import map_coordinates

class padfplot_dims():

    def __init__(self,rmin=0.0,rmax=1.0,rval=0.5,rval2=0.5,thval=0.0,rwid=0.05,thmin=0,thmax=360):
        
        self.rmin = rmin
        self.rmax = rmax
        self.rval = rval
        self.rval2 = rval2
        self.thval = thval
        self.rwid = rwid
        self.thmin = thmin
        self.thmax = thmax

    def get_ir(self, nr, r, flt=False):
        ir = nr*(r-self.rmin)/(self.rmax-self.rmin)
        if flt==False: ir = int(ir)
        return ir

    def get_ith(self, nth, th, flt=False):
        ith = nth*(th-self.thmin)/(self.thmax-self.thmin)
        if flt==False: ith=int(ith)
        return ith

    def get_rwid_index(self,nr,flt=False):
        irwid = int(nr*(self.rwid)/(self.rmax-self.rmin))
        if flt==False: irwid=int(irwid)
        return irwid    

# gamma scale

# multiply by power of r**n
def mult_radial_polynomial( corr, n, rmin=0, rmax=1 ):
    """multiply the correlation volume by r**n    
    """
    r = (rmax-rmin)*np.arange(corr.shape[0])/corr.shape[0] + rmin
    
    for i in np.arange(corr.shape[0]):
        for j in np.arange(corr.shape[0]):
            if (i !=0)and(j!=0):
                 corr[i,j,:] *= (r[i]*r[j])**n
    return corr


# remove angular average
def remove_angular_average( corr ):
    """subtract the angular average at each value of r and r prime (q and q prime)
    """
    for i in np.arange(corr.shape[0]):
        for j in np.arange(corr.shape[0]):
                 corr[i,j,:] += -np.average(corr[i,j,:])
    return corr


# multiply by |sin\theta|
def multiply_by_sintheta( corr, thmin, thmax):

    th = (np.pi/180.0)*(thmax-thmin)*np.arange(corr.shape[2])/corr.shape[2]
    absth = np.abs(np.sin(th))
    for i in np.arange(corr.shape[0]):
        for j in np.arange(corr.shape[0]):
                 corr[i,j,:] *= absth
    return corr

    


# generate unit labels (replacing Angstrom symbols)
def generate_unit_labels( runits,rq,rscale ):


    if runits in ['A','Angstrom']:
            ustr = u'\u00c5'
    elif runits in  ['um', 'micron']:
            ustr = r'$\mu$m'
    else:
            ustr = runits

    if rq=='q':
            ustr += r'$^{-1}$'    

    if rscale!=1:
        ustr = " x "+str(rscale)+" "+ustr

    rlabel = rq+r' ('+ustr+')'
 
    thlabel = r'$\theta$ (degrees)'

    return rlabel, thlabel


# make a 3D Gaussian
def make_gaussian_n(nx, ny, nz, rad=None, rady=-1., radz=-1., cenx=None, ceny=None, cenz = None, invert=0, norm=False ): 
    #set defaults
    if rad is None: rad = np.min(nx,ny)/2
    if cenx is None: cenx = nx//2
    if ceny is None: ceny = ny//2
    if cenz is None: cenz = nz//2
    radsq = rad**2
    if rady == -1.:
        radysq = radsq
    else:
        radysq = rady**2

    if radz == -1.:
        radzsq = radsq
    else:
        radzsq = radz**2
        
    xyz = np.mgrid[ -nx//2:nx//2-1:nx*1j, -ny//2:ny//2-1:ny*1j, -nz//2:nz//2-1:nz*1j]
    x, y, z = xyz[0], xyz[1], xyz[2]
    
    a = np.zeros([nx,ny,nz])
    a = np.exp(-(x**2)/radsq  - ( y**2)/radysq - ( z**2)/radzsq)
    #print( a.shape, x.shape, nx )
    a[ nx//2, ny//2, nz//2 ] = 1.0

    a = np.roll( a, cenx-nx//2, 0 )
    a = np.roll( a, ceny-ny//2, 1 )
    a = np.roll( a, cenz-nz//2, 2 )

    # normalise if required
    if norm == True: a *= 1./np.sum(a)
    
    return a
#end make_gaussian

# convolve by a 3D gaussian
def convolve_gaussian_n( image, rad=3, rady=1, radz=1):
     print( image.shape )
     c = make_gaussian_n( image.shape[0], image.shape[1], image.shape[2], rad, rady, radz, cenx=0, ceny=0, cenz=0 )
     fc = np.fft.fftn( c )
     fimage = np.fft.fftn( image )
     output = np.fft.ifftn( np.conjugate(fc)*fimage )
     return output

# scale the correlation slice by gamma
def corr_rescale( plane, gamma ):

     disp = plane*0.0
     ihigh = np.where( plane > 0 )
     ilow = np.where( plane < 0 )
     disp[ihigh] = np.abs( plane[ihigh] )**gamma
     disp[ilow] = - np.abs( plane[ilow] )**gamma

   #  disp[ihigh] = np.log(np.abs( plane[ihigh] ))
   #  disp[ilow] = - np.log(np.abs( plane[ilow] ))

     return disp

def costh_coordinate_correction(disp,ppdims):
     
    dispout = disp*0.0

    if disp.ndim==1:
        nth = disp.shape[0]
        th = ppdims.thmin + (ppdims.thmax-ppdims.thmin)*np.arange(nth)/nth
        z = nth*(np.cos(th[::-1]*np.pi/180)+1)/2            
        dispout = map_coordinates(disp, [z], order=3)
    elif disp.ndim==2:
        nth = disp.shape[1]
        th = ppdims.thmin + (ppdims.thmax-ppdims.thmin)*np.arange(nth)/nth
        z = nth*(np.cos(th[::-1]*np.pi/180)+1)/2
        for i in np.arange( disp.shape[0] ):
            zcorr = map_coordinates(disp[i,:], [z], order=3)
            dispout[i,:] = zcorr
    return dispout


# extract slice (see below for options)
def extract_section( corr, ppdims, stype='reqr', csampling=False ):
    #stype: 'reqr', 'rconst', 'roffset', 'thconst'
    print("Section extracted :",stype)

    if stype=='reqr':     
        disp = np.zeros( (corr.shape[0], corr.shape[2] ) )
        for i in np.arange(corr.shape[0]):
             disp[i,:] = corr[i,i,:]
        if csampling: disp = costh_coordinate_correction(disp,ppdims)

    elif stype=='rconst': 
        ir = ppdims.get_ir(corr.shape[0],ppdims.rval)
        disp = np.real(corr[:,ir,:])
        if csampling: disp = costh_coordinate_correction(disp,ppdims)

    elif stype=='thconst':        
        ith = ppdims.get_ith(corr.shape[2],ppdims.thval)
        disp = np.real(corr[:,:,ith])

    elif stype=='thline': 
        ir  = ppdims.get_ir(corr.shape[0],ppdims.rval)
        ir2 = ppdims.get_ir(corr.shape[0],ppdims.rval2)
        w = ppdims.get_rwid_index(corr.shape[0])
        disp = np.real(np.sum(np.sum(corr[ir-w:ir+1+w,ir2-w:ir2+1+w,:],0),0))
        if csampling: disp = costh_coordinate_correction(disp,ppdims)

    elif stype=='rline':        
        ith = ppdims.get_ith(corr.shape[2],ppdims.thval)
        disp = np.zeros(corr.shape[0])
        for i in np.arange(corr.shape[0]):
            disp[i] = corr[i,i,ith]

    return disp
    

