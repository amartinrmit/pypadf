
"""
paramsDIFFCORR.py

set up a params class to read/store the input parameters for a diffract_and_correlate calculation

Classes:
    paramsDIFFCORR - params class instance with diffraction and correlation input parameters
"""

from params.params import params
import os
import glob
import numpy as np

#
# class to set up the PADF parameters
#


class paramsDIFFCORR(params):
    """
    params class instance with PADF input parameters  

    Input parameter names - see help() for information about a specific parameter:
    config, outpath, tag, samplepath, nx, ny, nl, nlmin, nr, nq, nth, qmin, qmax, rmax, wl
    nthreads, dz, pw, nstart, npatterns, maskflag, maskname, rbin, diffcorrflag, outputdpflag,
    dpshiftflag, nxcrop, nycrop, shiftx, shifty, bg_estimate
    Parameters
    ----------
    config : str
        name of the config file with parameter values defined

    samplepath : str
        directory of diffraction files and file name format
        e.g. "/data/mydirectory/*.npy"
        file formats supported: .npy, .tiff, .dbin

    outpath : str
        path to save output of the padf calculation

    tag : str
        prefix for all output file names
        (used not to overwrite previous calculations)

    method : str
        method for matrix inversion 'svd' or 'legendre'
        'svd' inverts matrix with singular value decomposition
        'legendre' uses orthogonality properties of legendre polynomials

    nx : int
`       number of rows in input diffraction array

    ny : int
        number of columns in input diffraction array

    nl : int
        maximum order of spherical harmonics
        This sets the angular resolution of the padf

    nlmin : int
        minimum order of spherical harmonics to use

    nr : int
        number of radial samples in real space

    nq : int
        number of radial samples in q-space (must match input correlation file)
        this is essential when using dbin format

    nth : int
        number of angular samples
        must be the same in q-space and real-space

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

    nthreads : int
        number of CPU processes to use

    dz : float
        sample to detector distance. Units: metres

    pw : float
        width of detector pixel. Units: metres

    nstart : int
        index of first diffraction pattern to process (default 0)
    
    npatterns: int
        number of diffraction patterns to process

    maskflag : bool
        Correct correlation with mask correlation function (True/False)
    
    maskname : str
        path and filename of mask image to use (if maskflag is True)
        Mask takes values 0 for pixels to ignore and 1 for pixels to include.

    rbin : int
        rebinning factor. Must be a power of 2 (e.g. 2, 4, 8, 16)
        if set to 1, then diffraction patterns are not rebinned

    diffcorrflag : bool
        Calculate a difference correlation (True/False)    

    outputdpflag : bool
        Write out diffraction patterns to file (True/False)

    dpshiftflag : bool
        Shift the centre of the diffraction patterns (True/False)

    shiftx : float
        amount (in pixels) to shift centre of diffraction pattern (row coordinate)
        Only used if dpshiftflag is True.
        Non-integer shifts are interpolated.
    
    shifty : float
        amount (in pixels) to shift centre of diffraction pattern (column coordinate)
        Only used if dpshiftflag is True.
        Non-integer shifts are interpolated.

    nxcrop : int
        number of rows in cropped diffraction pattern.
        Diffraction pattern will be cropped if nxcrop and nycrop are positive integers.

    nycrop : int
        number of columns in cropped diffraction pattern.
    
    bgestimate : bool
        Calculate the background correlation (True/False)
    
    npulse : float
        estimated number of photons per exposure (per pulse)
        
    beamarea : float
        estimated beam area

    poisson : bool
        if true, the diffraction pattern will be sampled with a Poisson distribution, generating photon statistics

    polarisation : bool
        if true, diffraction pattern will be multiplied by a polarisation factor

    polx : float
        component of the polarisation vector in the horizontal direction

    poly : float
        compnonent of the polarisation vector is the vertical direction
    """


    def __init__(self):
        """
        Constructs an instance of the paramsCORR class. 
        Parameters are added using the add_parameter() method (derived from the params class)
        """
        
        # Correlation and common parameters
        ch = ["PY3CORRELATION"]

        params.__init__(self,"Correlation parameters",configheaders=ch)


        self.add_parameter("config", "None", cmdline="--config",cmdline2="-c", 
                           help="Name of file with all input parameters",
                           nargs=1,header=ch[0],pathflag=True)


        self.add_parameter("outpath", "None", cmdline="--outpath",cmdline2="-o",
                           help="Path where files will be written.",
                           nargs=1,header=ch[0],pathflag=True)
        
        self.add_parameter("tag", "None", cmdline="--tag",cmdline2="-t", 
                           help="text string prefix for each output file.",
                           nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("samplepath", "None", cmdline="--samplepath",cmdline2="-sp", 
                   help="Path or filelist where diffraction files are located.",
                   nargs=1,header=ch[0],pathflag=True)

        self.add_parameter("nq", int(100), cmdline="--nq",cmdline2="-nq", 
                           help="Number of q-space radial samples",
                           nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("nx", int(-1), cmdline="--nx",cmdline2="-nx", 
                           help="Number of x samples",
                           nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("ny", int(-1), cmdline="--ny",cmdline2="-ny", 
                           help="Number of y samples",
                           nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("nth", int(100), cmdline="--nth",cmdline2="-nth", 
                           help="Number of theta samples",
                           nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("nthreads", int(100), cmdline="--nthreads",cmdline2="-n", 
                           help="Number of processes to use",
                           nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("nstart", int(0), cmdline="--nstart",cmdline2="-ns", 
                           help="Number of starting diffraction pattern",
                           nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("npatterns", int(100), cmdline="--npatterns",cmdline2="-np", 
                           help="Number of patterns to correlate",
                           nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("rebin", int(0), cmdline="--rebin",cmdline2="-r", 
                           help="Integer factor to rebin by (multiple of 2)",
                           nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("wl", 0.0, cmdline="--wl",cmdline2="-w", help="wavelength (A)",
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("dz", 0.0, cmdline="--dz", cmdline2="-z",help="detector distance (m)",
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("pw", 0.0, cmdline="--pw", cmdline2="-pw",help="pixel width (m)",
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("bg_estimate", False, cmdline="--bg_estimate",cmdline2="-bg", 
                        help="flag to set correlation background calculation",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("diffcorrflag", False, cmdline="--diffcorrflag",cmdline2="-dc", 
                 help="flag to set difference correlation calculation",
                 nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("outputdp", False, cmdline="--outputdp",cmdline2="-od",
                 help="flag to write diffraction patterns to file",
                 nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("maskflag", False, cmdline="--maskflag",cmdline2="-mf", 
                 help="flag to peform mask correlation",
                 nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("maskname", "None", cmdline="--maskname",cmdline2="-mk",
                           help="File path to mask file.",
                           nargs=1,header=ch[0],pathflag=True)
        
        self.add_parameter("dp_shift_flag", False, cmdline="--dp_shift_flag",cmdline2="-dsf", 
                 help="set this flag to crop diffraction patterns",
                 nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("cropflag", False, cmdline="--cropflag",cmdline2="-cf",
                 help="set this flag to crop diffraction patterns",
                 nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("nxcrop", -1, cmdline="--nxcrop",cmdline2="-xc",
                 help="no. x pixels after cropping",
                 nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("nycrop", -1, cmdline="--nycrop",cmdline2="-yc", 
                 help="no. y pixels after cropping",
                 nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("shiftx", 0.0, cmdline="--shiftx", cmdline2="-sx",
                 help="shift diffraction pattern in x direction (no. of pixels)",
                 nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("shifty", 0.0, cmdline="--shifty",cmdline2="-sy", 
          help="shift diffraction pattern in y direction (no. of pixels)",
          nargs=1,header=ch[0],pathflag=False)


        self.add_parameter("chunksize", -1, cmdline="--chunksize",cmdline2="-cs", 
                 help="number of diffraction patterns to calculate at a time",
                 nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("writefreq", -1, cmdline="--writefreq",cmdline2="-wf",
                 help="write out correlation function at this freqency",
                 nargs=1,header=ch[0],pathflag=False)

        # diffraction parameters 
        self.add_parameter("pdbname", "None", cmdline="--pdbname",cmdline2="-po", help="Name of the pdb file containing atomic coordinates",
                nargs=1,header=ch[0],pathflag=True)
        self.add_parameter("diffoutpath", "None", cmdline="--diffoutpath",cmdline2="-do", help="Path where diffraction files will be written.",
                        nargs=1,header=ch[0],pathflag=True)
        self.add_parameter("fext", "None", cmdline="--fext",cmdline2="-x", help="extension for 2D diffraction output.",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("cenx", int(-1), cmdline="--cenx",cmdline2="-cx", help="x-coordinate (pixel units) of beam center)",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("ceny", int(-1), cmdline="--ceny",cmdline2="-cy", help="y-coordinates (pixel units) of beam center)",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("henkeflag", False, cmdline="--henkeflag",cmdline2="-hf", help="Boolean to set use of Henke parameters",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("rx", 1.0, cmdline="--rx",cmdline2="-rx", help="x component of rotation axis (arb units)",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("ry", 0.0, cmdline="--ry",cmdline2="-ry", help="y component of rotation axis (arb units)",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("rz", 0.0, cmdline="--rz",cmdline2="-rz", help="z component of rotation axis (arb units)",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("rtheta", 0.0, cmdline="--rtheta",cmdline2="-th", help="rotation angle (in degrees)",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("rotflag", False, cmdline="--rotflag",cmdline2="-rf", help="Rotate sample or not - Boolean",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("rmax", -1.0, cmdline="--rmax",cmdline2="-rmax", help="crop atoms within circular volume",
                        nargs=1,header=ch[0],pathflag=False) 

        self.add_parameter("alen", -1.0, cmdline="--alen",cmdline2="-al", help="width of the sample volume - x",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("blen", -1.0, cmdline="--blen",cmdline2="-bl", help="width of the sample volume - y",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("clen", -1.0, cmdline="--clen",cmdline2="-cl", help="width of the sample volume - z",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("periodicshift", False, cmdline="--periodicshift",cmdline2="-psh", help="Random periodic translation, \
                        using alen/blen/clen values - Boolean",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("npulse", 1e11, cmdline="--npulse",cmdline2="-nph", help="number of photons per pulse (exposure)",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("beamarea", 1e-14, cmdline="--beamarea",cmdline2="-ba", help="bream area in m^2",
                        nargs=1,header=ch[0],pathflag=False) 

        self.add_parameter("poisson", False, cmdline="--poisson",cmdline2="-ps", help="If true, diffraction pattern is sampled with the Poisson distribution to produce photon counts",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("polarisation", False, cmdline="--polarisation",cmdline2="-pol", help="bream area in m^2",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("polx", 1, cmdline="--polx",cmdline2="-px", help="component of polarisation vector in the horizontal direction",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("poly", 0, cmdline="--poly",cmdline2="-py", help="component of polarisation vector in the vertical direction",
                        nargs=1,header=ch[0],pathflag=False)

 

    def read_config_file(self):
        """Read the values of the input parameters from a text (config) file.
        
        All parameters are written to a log file in the outpath.
        """

        #cwd = os.getcwd()+"/"
        #self.read_parameters_from_file(self.parser.parse_args().config[0] )
        abspath = os.path.abspath(self.parser.parse_args().config[0] )
        self.parse_config_file(abspath)
        print("config file name:", abspath )
        self.qmax_calc()
        self.add_parameter( "qmax", self.qmax, cmdline="--dontuse", 
                            help="dont input qmax is calculated by the code", 
                            nargs=1, header="DONTUSE", pathflag=False )
        outname = self.makefname( self.outpath, self.tag, "_diffcorr_parameter_log", ".txt")
        self.write_params_to_file( outname )


        rnorm = np.sqrt( self.rx*self.rx + self.ry*self.ry + self.rz*self.rz)
        if rnorm > 0.0:
            self.rx *= 1.0/rnorm
            self.ry *= 1.0/rnorm
            self.rz *= 1.0/rnorm
        self.rtheta *= np.pi/180.0

    def load_flist_from_samplepath(self):
        """List the diffraction files in the samplepathi
            
           List is stored in self.flist
        """

        
        base = os.path.splitext(self.samplepath)
        
        if len(base[1])==0:
            self.flist = glob.glob(base[0]+"*.*")
        elif (base[1]==".txt") or (base[1]==".list"):
            self.flist = []
            for line in open(self.samplepath,'r'):
                self.flist.append(line)
        else:
            try:
                self.flist = glob.glob(self.samplepath)
            except:
                print("Check the sample path - not a valid directory, \
                      filelist or regular expression")

    def qmax_calc(self):
        """Calculate the maximum q value at the edge of the detector (or diffraction pattern)
           
            Assumes diffraction pattern is centred in the array
        """
        thmax = np.arctan( (self.nq)*self.pw/self.dz)
        self.qmax = (2/self.wl) * np.sin( thmax/2.0 )
