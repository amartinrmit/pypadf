
"""
paramsMASK.py

set up a params class to read/store the input parameters for masking a correlation function

Classes:
    paramsMASK - params class instance with maskcorr.py input parameters
"""

from params.params import params
import os


class paramsMASK(params):
    """
    params derived class with maskcorr input parameters  

    Parameters
    ----------
    config : str
        name of the config file with parameter values defined

    corrfile : str
        name of the q-space intensity correlation function
        supported formats : dbin, npy

    outpath : str
        path to save output of the padf calculation

    suffix : str
        suffix for all output file names
        (used not to overwrite previous calculations)

    submean : Bool
        subtract the angular mean from each q/q' line of the correaltion function (Default=True)
 
    sintheta : Bool
        mutliply the correlation function by abs(sin\theta) (Default=True)

    maskq : Bool
        set correlation function at high and/or low q values to zero (Default=False)

    maskth : Bool
        set correlation function near theta=0 to zero (Default=False)

    qmin : float
        minimum value of q in the correlation file. Units: inverse metres 
        (NOTE q IS DEFINED WITHOUT A FACTOR OF 2PI CONTRARY TO USUAL X-RAY CONVENTIONS)

    qmax : float
        maximum value of q in the correlation file. Units: inverse metres 
        (NOTE q IS DEFINED WITHOUT A FACTOR OF 2PI CONTRARY TO USUAL X-RAY CONVENTIONS)
    
    qmaskhigh : float
        upper limit of the q mask. Default value = qmax. Units: inverse metres 

    qmasklow : float
        lower limit of the q mask. Default value = qmin. Units: inverse metres

    thlim : float
        half-width of theta range to mask centred on theta=0. Units: degrees     
    
    thlimnorm : float
        bounds range for calculating the theta average. Ignored if less than 0. Units: degrees
    
    """


    def __init__(self):
        """Construct the paramsMASK class instance. 
           Adds all parameters to the class with default values.
        """        


        ch = ["MASK"]

        params.__init__(self,"Mask parameters",configheaders=ch)


        self.add_parameter("config", "None", cmdline="--config",cmdline2="-c", help="Name of file with all input parameters",
                nargs=1,header=["PADF"],pathflag=True)

        self.add_parameter("corrfile", "None", cmdline="--corrfile",cmdline2="-cf", help="Input q-space Correlation file.",
                        nargs=1,header=ch[0],pathflag=True)

        self.add_parameter("outpath", "None", cmdline="--outpath",cmdline2="-o", help="Path where files will be written.",
                        nargs=1,header=ch[0],pathflag=True)
        
        self.add_parameter("suffix", "processed", cmdline="--suffix",cmdline2="-s", help="text string suffix for each output file.",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("submean", True, cmdline="--submean",cmdline2="-sm", help="Set to subtract angular mean from each q/q' ring in the correlation function",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("sintheta", True, cmdline="--sintheta",cmdline2="-sth", help="Scale the correlation function by |sin\theta|",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("maskq", False, cmdline="--maskq",cmdline2="-mq", help="Set to mask q regions (set to zero)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("maskth", False, cmdline="--maskth",cmdline2="-mth", help="Set to mask angular regions (set to zero)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("qmin", 0.0, cmdline="--qmin",cmdline2="-qn", help="minimum q value (A^-1)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("qmax", 10.0, cmdline="--qmax",cmdline2="-qx", help="maximum q value (A^-1)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("qmaskhigh", 10.0, cmdline="--qmaskhigh",cmdline2="-qmh", help="upper bound on q mask (A^-1)",        
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("qmasklow", 10.0, cmdline="--qmasklow",cmdline2="-qml", help="lower bound on q mask (A^-1)",        
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("thlim", 10, cmdline="--thlim",cmdline2="-thx", help="bound on theta mask near zero (radians)",        
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("thlimnorm", -1.0, cmdline="--thlimnorm",cmdline2="-thxnorm", help="bound on theta range near zero for subtracting mean (radians)",        
                        nargs=1,header=ch[0],pathflag=False)


    def read_config_file(self):
        """Read the values of the input parameters from a text (config) file.
        
        All parameters are written to a log file in the outpath.
        """

        #cwd = os.getcwd()+"/"
        #self.read_parameters_from_file(self.parser.parse_args().config[0] )
        #print("config file name:", self.parser.parse_args().config[0] )
        abspath = os.path.abspath(self.parser.parse_args().config[0])
        self.parse_config_file( abspath )
        print("config file name:", abspath )
        #outpath = self.path_to_string(self.outpath)
        cf = self.path_to_string(self.corrfile)
        self.corrfile_basenoext = os.path.basename(os.path.splitext(cf)[0])
        logname = self.makefname( self.outpath, self.corrfile_basenoext, self.suffix+"_maskcorr_parameter_log", ".txt")
        self.write_params_to_file( logname )


#    def makefname( self, path, tag, suffix, fext):
#            outname = str(path.resolve())+tag+suffix+fext
#            #outname = path+tag+suffix+fext
#            return outname
