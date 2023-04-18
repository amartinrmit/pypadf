
"""
paramsPLOT.py

set up a params class to read/store the input parameters for plotting a correlation or a padf function

Classes:
    paramsPLOT - params class instance with plotfxs3D.py input parameters
"""

from params.params import params
import os
import numpy as np

class paramsPLOT(params):
    """
    params derived class with plotfxs3D.py input parameters  

    Parameters
    ----------
    config : str
        name of the config file with parameter values defined

    fname : str
        name of the q-space or real-space volume to display
        supported formats : dbin, npy

    outpath : str
        path to save output plots and figures as images

    suffix : str
        suffix for all output file names
        (used not to overwrite previous calculations)

    rmin : float
        minimum value of r (or q) in the input volume

    rmax : float
        maximum value of r (or q) in the input volume

    rq : string
        set real or qspace (affects axis labels and units)

    runits : string
        radial units. string used as entered, except 'Angstrom' will be replaced by symbol.
    
    thmin : float 
        minimum value of theta in the input volume (degrees)

    thmax : float
        maximum value of theta in the input volume (degrees)

    power : float
        multiply the radial dimension by r^power (or q^power)

    sintheta : Bool
        multiply volume by |sin\theta| before displaying images

    submean : Bool
        subtract the angular mean from each r/r' line

    convolve : Bool
        convolve the volume with a gaussian

    rwid : float
        radial width for gaussian convolution

    thwid : float
        angular width for gaussian convolution

    stype : str
        type of seciton to plot: reqr, rconst, thconst, rline, thline        

    rval : float
       value of r used in plotting (plotconstr / plotroffset / plotangleline)

    rpval : float
       value of r prime (plotangleline)

    thval : float
       value of theta (plotangleline / plotconsttheta)
        
    scale : float
       upper clim val in np.max(image)*scale

    scalel : float
       lower clim val is np.min(image)*scalel

    log : Bool
       take log of the images / line plots

    gamma : Bool
       gamma value to apply to images

    climhigh : float
        upper clim  on an absolute scale (takes priority over scale)

    climlow : float
        lower clim on an absolute scale (takes priority over scalel)

    """


    def __init__(self):
        """Construct the paramsPLOT class instance. 
           Adds all parameters to the class with default values.
        """        


        ch = ["FXSPLOT"]

        params.__init__(self,"Mask parameters",configheaders=ch)


        self.add_parameter("config", "None", cmdline="--config",cmdline2="-c", help="Name of file with all input parameters",
                nargs=1,header=["PADF"],pathflag=True)

        self.add_parameter("fname", "None", cmdline="--fname",cmdline2="-f",help="Input q-space or real-pace volume (r,r',\theta) or (q,q',theta).",
                        nargs=1,header=ch[0],pathflag=True)

        self.add_parameter("outpath", "None", cmdline="--outpath",cmdline2="-o", help="Path where files will be written.",
                        nargs=1,header=ch[0],pathflag=True)
        
        self.add_parameter("suffix", "processed", cmdline="--suffix",cmdline2="-s", help="text string suffix for each output file.",
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("submean", True, cmdline="--submean",cmdline2="-sm", help="Set to subtract angular mean from each q/q' ring in the correlation function",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("sintheta", True, cmdline="--sintheta",cmdline2="-sth", help="Scale the correlation function by |sin\theta|",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("convolve", True, cmdline="--convolve",cmdline2="-cv", help="Convolve the volume by a Gaussian",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("rmin", 0.0, cmdline="--rmin",cmdline2="-rn",
                         help="minimum radial value (r or q) in the input volume",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("rmax", 10.0, cmdline="--rmax",cmdline2="-rx", 
                        help="maximum radial value (r or q) in the input volume",
                        nargs=1,header=ch[0],pathflag=False)
        
        self.add_parameter("rq", 'r', cmdline="--rq", cmdline2="-rq",
                        help="set real-space 'r' or reciprocal space 'q'",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("rmindisp", 0.0, cmdline="--rmindisp",cmdline2="-rnd", 
                        help="minimum radial value (r or q) to display",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("rmaxdisp", 10.0, cmdline="--rmaxdisp",cmdline2="-rnx", 
                        help="maximum radial value (r or q) to display",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("runits", "m", cmdline="--runits",cmdline2="-ru", 
                        help="radial units to use on axis labels",        
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("thmin", 0.0, cmdline="--thmin", cmdline2="-thn",
                        help="minimum theta value in volume",        
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("thmax", 2*np.pi, cmdline="--thmax",cmdline2="-thx", 
                        help="maximum theta value in the volume",        
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("power", 0.0, cmdline="--power",cmdline2="-p", 
                        help="multiply the radial dimension by r**power or q**power",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("stype", 'reqr', cmdline="--stype",cmdline2="-y", 
                        help="type of section to plot: reqr, rconst, thconst, rline, thline",        
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("rval", 1.0, cmdline="--rval", cmdline2="-r",
                        help="r value to use for constant r slice or theta line",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("rpval", 1.0, cmdline="--rpval",cmdline2="-rp", 
                        help="r-prime value to use for theta line",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("thval", 0.0, cmdline="--thval",cmdline2="-th", 
                        help="theta value to use for constant theta slice or radial line (rline)",
                        nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("scale", 1.0, cmdline="--scale",cmdline2="-sc",help="upper clim val in np.max(image)*scale", nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("scalel", 1.0, cmdline="--scalel",cmdline2="-scl", help="lower clim val in np.min(image)*scalel", nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("clow", -1.0, cmdline="--clow",cmdline2="-cl", help="absolute clim lower limit (priotity over scale)", nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("chigh", -1.0, cmdline="--chigh",cmdline2="-ch", help="absolute clim upper limit (priority over scalel)", nargs=1,header=ch[0],pathflag=False)

        self.add_parameter("log", False, cmdline="--log",cmdline2="-l", help="log scale on all plots", nargs=1, header=ch[0],pathflag=False)

        self.add_parameter("gamma", -1.0, cmdline="--gamma",cmdline2="-g", help="gamma value to apply to all images", nargs=1, header=ch[0],pathflag=False)


        self.add_parameter("rwid", 1.0, cmdline="--rwid",cmdline2="-rw", help="radial width for convolution", nargs=1, header=ch[0],pathflag=False)
        self.add_parameter("thwid", 1.0, cmdline="--thwid",cmdline2="-thw", help="angular witdth for convolution", nargs=1, header=ch[0],pathflag=False)


    def read_config_file(self):
        """Read the values of the input parameters from a text (config) file.
        
        All parameters are written to a log file in the outpath.
        """
        #cwd = os.getcwd()+"/"
        #self.read_parameters_from_file(self.parser.parse_args().config[0] )
        #print("config file name:", self.parser.parse_args().config[0] )
        #self.fname_basenoext = os.path.basename(os.path.splitext(self.fname)[0])
        #self.write_params_to_file( self.outpath+self.fname_basenoext+"_"+self.stype+"_plot_parameter_log.txt" )

        abspath = os.path.abspath(self.parser.parse_args().config[0])
        self.parse_config_file( abspath )
        print("config file name:", abspath )
        fname = self.path_to_string(self.fname)
        self.fname_basenoext = os.path.basename(os.path.splitext(fname)[0])
        outname = self.makefname( self.outpath, self.fname_basenoext+"_"+self.stype, "_plot_parameter_log", ".txt")
        self.write_params_to_file( outname )

#    def makefname( path, tag, suffix, fext):
#            outname = str(path.resolve())+tag+suffix+fext
#            return outname
