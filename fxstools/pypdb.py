"""
pypdb.py

Read and process atomic information from a .pdb file

This script provides tools to read atomic and unit cell information from a 
pdb file. It extracts the information used by pydiffraction.py to calculate
diffraction patterns. It does not read all the pdb file information.

This file can not be executed independently. It is imported as a module 
and contains the following classes:

    * atomline  - represents properties of an atom
    * unitcell  - represents a unitcell
    * pdb       - read and store atom and unitcell information from a pdb file
"""

import numpy as np



class atomline:
    """
    A class to represent properties of an atom
    
    Attributes
    ----------
    elem : str
        1 or 2 letter symbol of the element
    atomN : int
        atomic number
    x : float
        x-coordinate
    y : float
        y-coordinate
    z : float
        z-coordinate
    occ : float
        occupancy
    Bfactor : float
        temperature factor
    """
    
    def __init__(self,elem="C",atomN=6, x=0.0,y=0.0,z=0.0,occ=1.0,Bfactor=1.0):
        """
        Constructs necessary attributes for an atom
        
        Parameters
        ----------
            elem : str
                1 or 2 letter symbol of the element
            atomN : int
                atomic number
            x : float
                x-coordinate
            y : float
                y-coordinate
            z : float
                z-coordinate
            occ : float
                occupancy
            Bfactor : float
                temperature factor
        """
        
        self.elem = elem
        self.atomN = atomN
        self.x = x
        self.y = y
        self.z = z
        self.occ = occ
        self.Bfactor = Bfactor


class unitcell:
    """
    A class to represent a unit cell
    
    Attributes
    ----------
    a : float
        length of lattice vector a
    b : float
        length of lattice vector b
    c : float
        length of lattice vector c
    alpha : float
        angle between b and c lattice vectors
    beta : float
        angle between c and a lattice vectors
    gamma : float
        angle between a and b lattice vectors
    volume: float
        volume of the unit cell
    """
    def __init__(self,a=1.0,b=1.0,c=1.0,alpha=90,beta=90,gamma=90):
        
        """
        Constructs attributes of a unit cell
        
        Attributes
        ----------
        a : float
            length of lattice vector a
        b : float
            length of lattice vector b
        c : float
            length of lattice vector c
        alpha : float
            angle between b and c lattice vectors
            input (degrees) auto converted to radians
        beta : float
            angle between c and a lattice vectors
            input (degrees) auto converted to radians
        gamma : float
            angle between a and b lattice vectors
            input (degrees) auto converted to radians
        volume: float
            volume of the unit cell
            
        Methods
        ---------
        deg2rad(theta)
            Converts theta into radians.
            
        volcalc()
            Calculates the volume of the unit cell.
            
        print_unit_cell()
            Prints the unit cell parameters to the console.
        
        """
        self.a = a
        self.b = b
        self.c = c
        self.alpha = self.deg2rad(alpha)
        self.beta = self.deg2rad(beta)
        self.gamma = self.deg2rad(gamma)
        self.volume = self.volcalc()

    def deg2rad( self, theta ):
        """Converts theta into radians.""" 
        return theta*np.pi/180.0

    def volcalc( self ):
        """Calculates the volume of the unit cell.

        Returns:
                 v (float) : volume of the unit cell
        """
        v = np.sqrt( 1 + 2*np.cos(self.alpha)*np.cos(self.beta)*np.cos(self.gamma)
                     - np.cos(self.alpha)**2 - np.cos(self.beta)**2 - np.cos(self.gamma)**2 )

        v *= self.a*self.b*self.c
        return v

    def print_unit_cell(self):
        """Prints the unit cell parameters to the console."""
        print( "a, b, c", self.a, self.b, self.c )
        print( "alpha, beta, gamma", self.alpha, self.beta, self.gamma )
        print( "unit cell volume", self.volume )

class pdb:
    """
    A class to read and store atom information from a pdb file
    Stores relevant information for pydiffraction calculations
     
    Attributes
    ----------
    pdbname : str
        name of the pdb file
    cryst : unitcell class object
        contains lattice parameters and volume
    atomlist : list of atomline objects
        stores atomic elements and positions
    all : dictionary of strs
        element codes and atomic numbers
    sing : dict of strs
        single letter element codes and atomic numbers
    dbl : dict of strs
        double letter element codes and atomic numbers
    xmin, xmax : float, float
        minimum and maximum x-coordiantes
    ymin, ymax : float, float
        minimum and maximum y-coordiantes
    zmin, zmax : float, float
        minimum and maximum z-coordiantes
    mvol : float
        rectangular volume defined by min/max coordinates 
        
    Methods
    ---------
    read_pdb()
        Reads a pdb file to extract atom and unit cell information
    
    set_up_element_codes()
        Defines dictionaries with element symbols and atomic numbers
    
    check_element_code(str)
        checks if str is one of the known element codes.
        
    sort_atom_list()
        Creates a nested list of atoms sorted according to element
        
    split_atom_line(line)
            splits an ATOM line of pdb file into strings for each parameter
            
    split_atom_line_no_strip( self, line):
        splits an ATOM line of pdb file into strings for each parameter
        without stripping white space.
        
    maxdims()
        Determines the minimum/maximum coordinate values from the atom list
        Calcualtes a volume from min/max coordinates.
    """
    def __init__(self,pdbname):
        """
        Constructor for pdb class
        Stores relevant atom and unitcell information for pydiffraction 
        calculations
         
        Attributes
        ----------
        pdbname : str
            name of the pdb file
        cryst : unitcell class object
            contains lattice parameters and volume
        atomlist : list of atomline objects
            stores atomic elements and positions
        all : dictionary of strs
            element codes and atomic numbers
        sing : dict of strs
            single letter element codes and atomic numbers
        dbl : dict of strs
            double letter element codes and atomic numbers
        xmin, xmax : float, float
            minimum and maximum x-coordiantes
        ymin, ymax : float, float
            minimum and maximum y-coordiantes
        zmin, zmax : float, float
            minimum and maximum z-coordiantes
        mvol : float
            rectangular volume defined by min/max coordinates 
        
        """
        self.pdbname = pdbname
        self.set_up_element_codes()
        self.cryst = None
        
        
    def read_pdb(self):
        """Reads atom and lattice information from pdb file
        
        Attibutes required
        ----------
        self.pdbname - name of the pdb file
        self.all     - dictionary of element codes and atomic numbers
        
        Attributes overwritten
        ----------
        The following attributes are defined or overwritten
        with information from the pdb file:
        
        self.atomlist  - list of atomic coordinates and elements
        self.cryst     - unit cell parameters
        """
        f = open(self.pdbname,"r")

        self.atomlist = []

        for line in f:

            bits = line.split()

            if bits[0] == "CRYST1":
                #print( line )
                self.cryst = unitcell( float(bits[1]), float(bits[2]), float(bits[3]),
                                  float(bits[4]), float(bits[5]), float(bits[6]))


            if bits[0] == "ATOM" or bits[0]=="HETATOM":

                b = self.split_atom_line( line )
                                
                elem = self.check_element_code( b[2] )
#                atom = atomline( elem, self.all[elem], 
#                                 float(bits[6]), float(bits[7]), float(bits[8]) )

                if elem in self.all.keys():
                    
                    atom = atomline( elem, self.all[elem], 
                                     float(b[7]), float(b[8]), float(b[9]),
                                     float(b[10]), float(b[11]) )

                    self.atomlist.append( atom )
                else:
                    print("Could not identify element code in pdb file:",elem)





        f.close()

        # sort the atom list to find the number of elements
        self.sort_atom_list()


    def write_pdb(self, outname):
    
        alstart = "ATOM  {num:5d}"
        xyzbit = "{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}"

        f = open(outname,'w')
        
        for i, al in enumerate(self.atomlist):
            if i!=0:f.write("\n")
            strng = alstart.format(num=i+1)+al.elem.rjust(3)+"XXX".rjust(6)+"X".rjust(2)+"{n:4d}    ".format(n=i+1)
            strng += xyzbit.format(x=al.x,y=al.y,z=al.z,occ=al.occ,b=al.Bfactor)+al.elem.rjust(12)
            f.write(strng)
        f.close()
        """
        if len(line)>=6: b.append( line[0:6] )      #line type 6 characters

        if len(line)>=11: b.append( line[6:11].strip() )     #atom number
        if len(line)>=16: b.append( line[13:16].strip() )    #element code
        if len(line)>=19: b.append( line[16:20].strip() )    # residue
        if len(line)>=21: b.append( line[20:22].strip() )
        if len(line)>=25: b.append( line[22:26].strip() )
        if len(line)>=29: b.append( line[26:30].strip() )
        if len(line)>=37: b.append( line[30:38].strip() )    # x
        if len(line)>=45: b.append( line[38:46].strip() )    # y
        if len(line)>=52: b.append( line[46:54].strip() )    # z
        if len(line)>=58: b.append( line[54:60].strip() )    # occupancy
        if len(line)>=64: b.append( line[60:66].strip() )    # B factor
        if len(line)>=77: b.append( line[66:79].strip() )    # final letter

        self.elem = elem
        self.atomN = atomN
        self.x = x
        self.y = y
        self.z = z
        self.occ = occ
        self.Bfactor = Bfactor
        """
                
    def set_up_element_codes(self):
        """Sets up dictionaries of element codes and atomic numbers
        
        Attibutes required
        ----------
        None
        
        Attributes overwritten
        ----------
        The following attributes are defined or overwritten
        with information from the pdb file:
        
        self.all  - dictionary of all atomic numbers and element codes
        self.sng  - dictionary of single letter element codes
        self.dbl  - dictionary of double letter element codes
        """

        self.all ={ "H":1,"He":2,"Li":3,"Be":4, "B":5, "C":6, "N":7, "O":8, "F":9,\
                   "Ne":10, "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16,
               "Cl":17,"Ar":18, "K":19, "Ca":20, "Sc":21, "Ti":22, "V":23, "Cr":24, "Mn":25,
               "Fe":26,"Co":27,"Ni":28,"Cu":29, "Zn":30,"Ga":31,"Ge":32,"As":33,"Se":34, 
               "Br":35,"Kr":36,"Rb":37, "Sr":38, "Y":39, "Zr":40, "Nb":41, "Mo":42, "Tc":43, "Ru":44,
               "Rh":45, "Pd":46, "Ag":47, "Cd":48, "In":49, "Sn":50, "Sb":51, "Te":52, "I":53, "Xe":54,
               "Cs":55,"Ba":56,
               "La":57,"Ce":58,"Pr":59,"Nd":60,"Pm":61,"Sm":62,"Eu":63,"Gd":64,"Tb":65,
               "Dy":66,"Ho":67,"Er":68,"Tm":69,"Yb":70,
               "Lu":71,"Hf":72,"Ta":73, "W":74,"Re":75,"Os":76,"Ir":77,
               "Pt":78,"Au":79,"Hg":80,"Tl":81,"Pb":82,"Bi":83,"Po":84,"At":85,"Rn":86,
               "Fr":87,"Ra":88, "Ac":89,"Th":90,"Pa":91, "U":92,"SP":0,"SH":0,"PP":0
               }

        self.sing = {}
        self.dbl = {}
        for k in self.all.keys():
            if len(k)==1:
                self.sing.update( {k : self.all[k]} )
            if len(k)==2:
                self.dbl.update( {k : self.all[k]} )
        

    def check_element_code(self,str):
        """Checks if a string matches an element code in the list
        
        Parameters
        ----------
        str - string to check against element code list
        
        Returns
        ----------
        Element code string, if the element code is found
        
        "NN-"+str, if the parameter str is not found in the element code list.
        """
        if str[:2] in self.dbl:
            return str[:2]
        elif str[0] in self.sing:
            return str[0]
        else:
            return "NN-"+str


    #
    # count the number of different types of elements in the pdb file
    #
    def sort_atom_list(self):
        """Creates a nested list of atoms sorted according to element
        
        Parameters
        ----------
        None
        
        Returns
        ----------
        None
        
        Attributes overwritten
        ----------
        self.elements - list of element codes found in the atom list
        self.zlist    - list of atomic numbers found in the atom list
        self.sorted_atom_list - nested list of sorted atoms. First index
            matches the order of self.elements
        """
        self.elements = []
        self.sorted_atom_list = []

        for a in self.atomlist:
            if a.elem not in self.elements:
                self.elements.append(a.elem)
                self.sorted_atom_list.append( [a] )
            else:
                ind = self.elements.index(a.elem)
                self.sorted_atom_list[ind].append(a)
                
        self.zlist = []
        for elem in self.elements:
            self.zlist.append( self.all[elem] )

    def split_atom_line( self, line):
        """splits an ATOM line of pdb file into strings for each parameter
        
        Parameters
        ----------
        line : str
            a string containing an ATOM line of the pdb file
        
        Returns
        ----------
        b : list of strs
            list of the parameters in the ATOM line, with whitespace stripped
        """
        b = []
        
        if len(line)>=6: b.append( line[0:6] )      #line type 6 characters

        if len(line)>=11: b.append( line[6:11].strip() )     #atom number
        if len(line)>=16: b.append( line[13:16].strip() )    #element code
        if len(line)>=19: b.append( line[16:20].strip() )    # residue
        if len(line)>=21: b.append( line[20:22].strip() )
        if len(line)>=25: b.append( line[22:26].strip() )
        if len(line)>=29: b.append( line[26:30].strip() )
        if len(line)>=37: b.append( line[30:38].strip() )    # x
        if len(line)>=45: b.append( line[38:46].strip() )    # y
        if len(line)>=52: b.append( line[46:54].strip() )    # z
        if len(line)>=58: b.append( line[54:60].strip() )    # occupancy
        if len(line)>=64: b.append( line[60:66].strip() )    # B factor
        if len(line)>=77: b.append( line[66:79].strip() )    # final letter

        return b


    def split_atom_line_no_strip( self, line):
        """splits an ATOM line of pdb file into strings for each parameter
           without stripping white space.
        
        Parameters
        ----------
        line : str
            a string containing an ATOM line of the pdb file
        
        Returns
        ----------
        b : list of strs
            list of the parameters in the ATOM line, with whitespace stripped
        """

        b = []
        
        if len(line)>=6: b.append( line[0:6] )      #line type 6 characters

        if len(line)>=11: b.append( line[6:11] )     #atom number
        if len(line)>=16: b.append( line[11:16] )    #element code
        if len(line)>=19: b.append( line[16:20] )    # residue
        if len(line)>=21: b.append( line[20:22] )
        if len(line)>=25: b.append( line[22:26] )
        if len(line)>=29: b.append( line[26:30] )
        if len(line)>=37: b.append( line[30:38] )    # x
        if len(line)>=45: b.append( line[38:46] )    # y
        if len(line)>=52: b.append( line[46:54] )    # z
        if len(line)>=58: b.append( line[54:60] )    # occupancy
        if len(line)>=64: b.append( line[60:66] )    # B factor
        if len(line)>=77: b.append( line[66:79] )    # final letter

        return b

    def maxdims( self ):
        """Determines the minimum/maximum coordinate values from the atom list
           Calcualtes a volume from min/max coordinates.
        
        Attibutes required
        ----------
        self.atomlist - list of atom information
        
        Attributes overwritten
        ----------
        The following attributes are defined or overwritten
        with information from the pdb file:
        
        self.xmin, self.xmax  - min/max x coordinate values
        self.ymin, self.ymax  - min/max y coordinate values
        self.zmin, self.zmax  - min/max z coordinate values
        self.mvol - volume of the rectangular volume defined by min/max values
        """
        xmin = 1e10
        xmax = -1e10
        ymin = 1e10
        ymax = -1e10
        zmin = 1e10
        zmax = -1e10
        

        for atom in self.atomlist:

            if atom.x < xmin:
                xmin = atom.x

            if atom.x > xmax:
                xmax = atom.x

            if atom.y < ymin:
                ymin = atom.y

            if atom.y > ymax:
                ymax = atom.y

            if atom.z < zmin:
                zmin = atom.z

            if atom.z > zmax:
                zmax = atom.z

        self.xmax, self.ymax, self.zmax = xmax, ymax, zmax
        self.xmin, self.ymin, self.zmin = xmin, ymin, zmin
        self.xlen, self.ylen, self.zlen = xmax-xmin, ymax-ymin, zmax-zmin                
        self.mvol = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
