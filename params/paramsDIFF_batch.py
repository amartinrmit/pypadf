
"""
paramsDIFF_batch.py

set up a params class to read/store the input parameters for a batch diffraction calculation

Classes:
    paramsDIFFBATCH - params class instance with batch diffraction input parameters
"""

from params.paramsDIFF import paramsDIFF
import os
import numpy as np

#
# class to set up the PADF parameters
#


class paramsDIFFBATCH(paramsDIFF):

    """
    params class instance with diffraction input parameters  

    Parameters
    ----------
    alen : float
        length of lattice parameter a
    
    blen : float
        length of lattice parameter b
   
    clen : float
        length of lattice parameter c
    """

    def __init__(self):

        ch = ["PYDIFFRACTION"]
        paramsDIFF.__init__(self)


        self.add_parameter("alen", 1.0, cmdline="--alen",cmdline2="-al", help="Length of lattice vector a, \
                        assumed cubic cell",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("blen", 1.0, cmdline="--blen",cmdline2="-bl", help="Length of lattice vector b, \
                        assumed cubic cell",
                        nargs=1,header=ch[0],pathflag=False)
        self.add_parameter("clen", 1.0, cmdline="--clen",cmdline2="-cl", help="Length of lattice vector c, \
                        assumed cubic cell",
                        nargs=1,header=ch[0],pathflag=False)



    def read_diff_parameters_from_file(self):
        """Read the values of the input parameters from a text (config) file.
        
        All parameters are written to a log file in the outpath.
        """

        #cwd = os.getcwd()+"/"
        #self.read_parameters_from_file(cwd+self.parser.parse_args().config[0] )
        abspath = os.path.abspath(self.parser.parse_args().config[0])
        self.parse_config_file( abspath )
        self.parse_commandline_args()
#        self.read_config_file(cwd+self.parser.parse_args().config[0] )
        print("config file name:", abspath )
        outname = self.makefname( self.outpath, self.tag, "_diffraction_batch_parameter_log", ".txt")
        self.write_params_to_file( outname )

        ## TODO: clean up the rotation axis - normalise vector etc; convert angle to radians etc
        rnorm = np.sqrt( self.rx*self.rx + self.ry*self.ry + self.rz*self.rz)
        if rnorm > 0.0:
            self.rx *= 1.0/rnorm
            self.ry *= 1.0/rnorm
            self.rz *= 1.0/rnorm
        self.rtheta *= np.pi/180.0
