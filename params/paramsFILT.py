
"""
paramsFILT.py

set up a params class to read/store the input parameters for a model_padf_filter.py calculation

Classes:
    paramsFILT - params class instance with filter input parameters
"""

from params.params import params
import os
import glob
import numpy as np

#
# class to set up the model-padf-filter parameters
#


class paramsFILT(params):
    """
    params class instance with model_padf_filter input parameters  

    Input parameter names - see help() for information about a specific parameter:
    config, outpath, tag, samplepath, nx, ny, nl, nlmin, nr, nq, nth, qmin, qmax, rmax, wl
    nthreads, dz, pw, nstart, npatterns, maskflag, maskname, rbin, diffcorrflag, outputdpflag,
    dpshiftflag, nxcrop, nycrop, shiftx, shifty, bg_estimate
    Parameters
    ----------
    config : str
        name of the config file with parameter values defined

    padfpath : str
        directory  and file name of input padf file

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

    """


    def __init__(self):
        """
        Constructs an instance of the paramsCORR class. 
        Parameters are added using the add_parameter() method (derived from the params class)
        """
        

        ch = ["PYBLFILTER"]

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

        self.add_parameter("padffile", "None", cmdline="--padffile",cmdline2="-pf", 
                   help="File name of the input padf file (with full path)",
                   nargs=1,header=ch[0],pathflag=True)

        self.add_parameter("nq", int(100), cmdline="--nq",cmdline2="-nq", 
                           help="Number of q-space radial samples",
                           nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("nr", int(-1), cmdline="--nr",cmdline2="-nr", 
                           help="Number of real space radial samples",
                           nargs=1,header=ch[0],pathflag=False)
        

        self.add_parameter("nth", int(100), cmdline="--nth",cmdline2="-nth", 
                           help="Number of theta samples",
                           nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("nl", int(32), cmdline="--nl",cmdline2="-nl",
                           help="Number of spherical harmonics",
                           nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("nlmin", int(0), cmdline="--nlmin",cmdline2="-nlm",
                           help="Minimum spherical harmonic l value",
                           nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("wl", 0.0, cmdline="--wl",cmdline2="-w", help="wavelength (m)",
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("qmin", 0.0, cmdline="--qmin",cmdline2="-qn", help="minimum q (m^-1)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("qmax", 0.0, cmdline="--qmax",cmdline2="-qx", help="maximum q (m^-1)",
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("rmax", 0.0, cmdline="--rmax",cmdline2="-rx", help="rmax (m)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("blqqtocorr", True, cmdline="--blqqtocorr",cmdline2="-bc", help="flag to map blqq to correlation file and back.",nargs=1,header=ch[0],pathflag=False) 

        self.add_parameter("interpolate", False, cmdline="--interpolate",cmdline2="-i", help="Use interpolation for q bins.",nargs=1,header=ch[0],pathflag=False) 

        self.add_parameter("order", 1, cmdline="--order",cmdline2="-or", 
                           help="order of interpolation (if used)",
                           nargs=1,header=ch[0],pathflag=False)



    def read_config_file(self):
        """Read the values of the input parameters from a text (config) file.
        
        All parameters are written to a log file in the outpath.
        """

        #cwd = os.getcwd()+"/"
        #self.read_parameters_from_file(cwd+self.parser.parse_args().config[0] )

        abspath = os.path.abspath(self.parser.parse_args().config[0])
        self.parse_config_file( abspath )
        print("config file name:", abspath )
        outpath = self.path_to_string(self.outpath)
        self.write_params_to_file( outpath+self.tag+"_blfilter_parameter_log.txt" )
