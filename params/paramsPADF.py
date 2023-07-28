
"""
paramsPADF.py

set up a params class to read/store the input parameters for a corrtopadf.py calculation

Classes:
    paramsPADF - params class instance with PADF input parameters
"""

from params.params import params
import os


class paramsPADF(params):
    """
    params derived class with PADF input parameters  

    Input parameter names - see help() for information about a specific parameter:
    config, corrfile, outpath, tag, method, nl, nlmin, nr, nq, nth, qmin, qmax,rmax, wl

    Parameters
    ----------
    config : str
        name of the config file with parameter values defined

    corrfile : str
        name of the q-space intensity correlation function
        supported formats : dbin, npy

    outpath : str
        path to save output of the padf calculation

    tag : str
        prefix for all output file names
        (used not to overwrite previous calculations)

    method : str
        method for matrix inversion 'svd' or 'legendre'
        'svd' inverts matrix with singular value decomposition
        'legendre' uses orthogonality properties of legendre polynomials

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

    density : float
        number density in atoms per Angstrom^-3 (default 0.1 A^-3). If not set, consider the padf in arbitraty units

    beamnorm: float
        normalization for the experimental intensity (r_e^2 I_0 d\Omega)
        (TODO add a script to pre-calculate)
    """


    def __init__(self):
        """Construct the paramsPADF class instance. 
           Adds all parameters to the class with default values.
        """        


        ch = ["PADF"]

        params.__init__(self,"PADF parameters",configheaders=ch)


        self.add_parameter("config", "None", cmdline="--config",cmdline2="-c", help="Name of file with all input parameters",
                nargs=1,header=["PADF"],pathflag=True)

        self.add_parameter("corrfile", "None", cmdline="--corrfile",cmdline2="-cf", help="Input q-space Correlation file.",
                        nargs=1,header=ch[0],pathflag=True)

        self.add_parameter("outpath", "None", cmdline="--outpath",cmdline2="-o", help="Path where files will be written.",
                        nargs=1,header=ch[0],pathflag=True)
        
        self.add_parameter("tag", "None", cmdline="--tag",cmdline2="-t", help="text string prefix for each output file.",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("method", "svd", cmdline="--method",cmdline2="-m", help="method to account for Ewald sphere: svd or legendre.",
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("nl", int(10), cmdline="--nl",cmdline2="-nl", help="Number of spherical harmonic values (l)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("nlmin", int(0), cmdline="--nlmin",cmdline2="-nlm", help="Minimum spherical harmonic value l",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("nr", int(100), cmdline="--nr",cmdline2="-nr", help="Number of real-space radial samples",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("nq", int(100), cmdline="--nq",cmdline2="-nq", help="Number of q-space radial samples",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("nth", int(100), cmdline="--nth",cmdline2="-nth", help="Number of theta samples",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("qmin", 0.0, cmdline="--qmin",cmdline2="-qn", help="minimum q value (A^-1)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("qmax", 10.0, cmdline="--qmax",cmdline2="-qx", help="maximum q value (A^-1)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("rmax", 10.0, cmdline="--rmax",cmdline2="-rx", help="maximum r value (A)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("wl", 0.0, cmdline="--wl",cmdline2="-w", help="wavelength (A)",
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("density", 0.1, cmdline="--density",cmdline2="-d", help="number density in (A^-3)",
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("beamnorm", 1.0, cmdline="--beamnorm",cmdline2="-b", help="dimensionless norm to convert intensity to structure factor squared",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("legendre_norm", True, cmdline="--lnorm",cmdline2="-ln", help="Normalise  the legendre polynomials",
                            nargs=1,header=ch[0],pathflag=False)   

    def read_config_file(self):
        """Read the values of the input parameters from a text (config) file.
        
        All parameters are written to a log file in the outpath.
        """

        #cwd = os.getcwd()+"/"
        #self.read_parameters_from_file(cwd+self.parser.parse_args().config[0] )
        abspath = os.path.abspath(self.parser.parse_args().config[0])
        self.parse_config_file( abspath )
        self.parse_commandline_args()
        self.checkpaths()
        print("config file name:", abspath )
        outname = self.makefname( self.outpath, self.tag, "_corrtopadf_parameter_log",".txt" )
        self.write_params_to_file( outname )
