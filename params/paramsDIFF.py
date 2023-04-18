
"""
paramsDIFF.py

set up a params class to read/store the input parameters for a diffraction calculation

Classes:
    paramsDIFF - params class instance with diffraction input parameters
"""

from params.params import params
import os
import numpy as np

#
# class to set up the PADF parameters
#


class paramsDIFF(params):
    """
    params class instance with diffraction input parameters  

    Parameters
    ----------
    config : str
        name of the config file with parameter values defined
    
    pdbname : str
        name of pdb file containing the atomic structure information

    outpath : str
        path to save output of the padf calculation

    tag : str
        prefix for all output file names
        (used not to overwrite previous calculations)
    
    fext : str
        file extension for the output diffraction pattern (e.g. .npy)    
 
    nx : int
`       number of rows in input diffraction array

    wl : float
        wavelength of the incident beam. Units: metres

    dz : float
        sample to detector distance. Units: metres

    pw : float
        width of detector pixel. Units: metres

    cenx : int
        row of the pixel location of the beam centre

    ceny : int
        column of the pixel location of the beam centre

    henkeflag : Bool
        Wavelength depedent correction to the atomic scattering factors (using Henke data)

    rx : float
        x component of the rotation axis (arb units).

    ry : float
        y component of the rotation axis (arb units).
    
    rz : float
        z component of the rotation axis (arb units)

    rtheta : float
        angle to rotate sample

    rotflag : Bool
        perform a sample rotation (True/False)

    npatterns : int
        number of diffraction patterns to calculate

    rmax : float
        crop atoms in a spherical volume centred on the origin (units to match pdb file)

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

        ch = ["PYDIFFRACTION"]
        params.__init__(self,"PYDIFFRACTION parameters",configheaders=ch)


        self.add_parameter("config", "None", cmdline="--config",cmdline2="-c", help="Name of file with all input parameters",
                nargs=1,header=ch[0],pathflag=True)

        self.add_parameter("pdbname", "None", cmdline="--pdbname",cmdline2="-p", help="Name of the pdb file containing atomic coordinates",
                nargs=1,header=ch[0],pathflag=True)
        self.add_parameter("outpath", "None", cmdline="--outpath",cmdline2="-o", help="Path where files will be written.",
                        nargs=1,header=ch[0],pathflag=True)
        self.add_parameter("tag", "None", cmdline="--tag",cmdline2="-t", help="text string prefix for each output file.",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("fext", "None", cmdline="--fext",cmdline2="-x", help="extension for 2D diffraction output.",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("nx", int(100), cmdline="--nx",cmdline2="-nx", help="Number of pixels on the side length of array",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("wl", 1e-10, cmdline="--wl",cmdline2="-w", help="Wavelength (metres)",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("dz", 0.1, cmdline="--dz",cmdline2="-z", help="Detector distance (metres)",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("pw", 0.00005, cmdline="--pw",cmdline2="-pw", help="Pixel width (metres)",
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
        self.add_parameter("rotflag", False, cmdline="--rotflag",cmdline2="-r", help="Rotate sample or not - Boolean",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("npatterns", 1, cmdline="--npatterns",cmdline2="-np", help="Number of patterns to simulate",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("rmax", -1.0, cmdline="--rmax",cmdline2="-rmax", help="crop atoms within circular volume",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("npulse", 1e11, cmdline="--npulse",cmdline2="-nph", help="number of photons per pulse (exposure)",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("beamarea", 1e-14, cmdline="--beamarea",cmdline2="-ba", help="beam area in m^2",
                        nargs=1,header=ch[0],pathflag=False) 

        self.add_parameter("poisson", False, cmdline="--poisson",cmdline2="-ps", help="If true, diffraction pattern is sampled with the Poisson distribution to produce photon counts",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("polarisation", False, cmdline="--polarisation",cmdline2="-pol", help="bream area in m^2",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("polx", 1, cmdline="--polx",cmdline2="-px", help="component of polarisation vector in the horizontal direction",
                        nargs=1,header=ch[0],pathflag=False) 
        self.add_parameter("poly", 0, cmdline="--poly",cmdline2="-py", help="component of polarisation vector in the vertical direction",
                        nargs=1,header=ch[0],pathflag=False)

 
    def read_diff_parameters_from_file(self):
        """Read the values of the input parameters from a text (config) file.
        
        All parameters are written to a log file in the outpath.
        """

        cwd = os.getcwd()+"/"
        print("config file name:", cwd+self.parser.parse_args().config[0],"\n" )
        #self.read_parameters_from_file(cwd+self.parser.parse_args().config[0] )
        self.parse_config_file(cwd+self.parser.parse_args().config[0] )
        #outname =  self.d["outpath"].value / (self.d["tag"].value+"_parameter_log.txt")
        outname = self.makefname( self.outpath, self.tag, "_diffract_parameter_log.txt")
        self.write_params_to_file( outname )

        rnorm = np.sqrt( self.rx*self.rx + self.ry*self.ry + self.rz*self.rz)
        if rnorm > 0.0:
            self.rx *= 1.0/rnorm
            self.ry *= 1.0/rnorm
            self.rz *= 1.0/rnorm
        self.rtheta *= np.pi/180.0
