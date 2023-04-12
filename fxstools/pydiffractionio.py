"""
pydiffractionio.py

Read and write image files and correlation volumes.

This file contains functions and can not be executed independently. 
It supports the following file types:

File reading: npy, bin, dbin, h5, tif
Image writing: dbin, npy

dbin = double precision binary file with no header information. 
       2D images are assumed to be square

bin = single precision binary file with no header information

"""
"""
Change log:

28/10/19 Jack:
Added tifread function using the python3 compatible PIL.Image
module for reading .tif images.

In read_image function, added explicit option for reading .tif files.

"""

import array
import imageio
import os
import struct
import h5py
from PIL import Image

#import amser as ser
import numpy as np


def read_dbin(fname, swapbyteorder=0, nx=-1, ny=-1):
    """Read an image from a double precision binary file withoutheader information
        
    Parameters
    ----------
    fname : string
        name of the double precision binary file
    
    swapbyteorder : integer
        1 (0) - do (not) swap byte order of each number
                may be required when reading/writing from different operating systems

    nx, ny : integer (optional)
        option to specify number of rows and columns in 2D image
 
    Returns
    ----------
    output : numpy array (float)
        the image or 2D data    
    """

    size = os.path.getsize(fname)
    # print "sizes:", size, size/8
    b = array.array('d')
    f = open(fname, "r")
    b.fromfile(f, size / 8)
    f.close()
    lst = b.tolist()

    n = len(lst)
    if (nx == -1) and (ny == -1):
        nx = int(round(n ** (1.0 / 2.0)))
    if ny == -1:
        ny = nx
    # print "nx", nx
    output = np.array(lst).reshape(nx, ny)
    if swapbyteorder == 1:
        output = output.newbyteorder()
    return output


def read_bin(fname, swapbyteorder=0, nx=-1, ny=-1):    
    """Read an image from a precision binary file withoutheader information
        
    Parameters
    ----------
    fname : string
        name of the double precision binary file
    
    swapbyteorder : integer
        1 (0) - do (not) swap byte order of each number
                may be required when reading/writing from different operating systems

    nx, ny : integer (optional)
        option to specify number of rows and columns in 2D image
 
    Returns
    ----------
    output : numpy array (float)
        the image or 2D data    
    """
    size = os.path.getsize(fname)
    # print "sizes:", size, size/8
    b = array.array('f')
    f = open(fname, "r")
    b.fromfile(f, size / 4)
    f.close()
    lst = b.tolist()

    n = len(lst)
    if (nx == -1) and (ny == -1):
        nx = int(round(n ** (1.0 / 2.0)))
    if ny == -1:
        ny = nx
    output = np.array(lst).reshape(nx, ny)
    if swapbyteorder == 1:
        output = output.newbyteorder()
    return output


def read_correlation(fname, swapbyteorder=0):
    """Read a double precision binary file without header information
      Does not reshape the data.
        
    Parameters
    ----------
    fname : string
        name of the double precision binary file
    
    swapbyteorder : integer
        1 (0) - do (not) swap byte order of each number
                may be required when reading/writing from different operating systems

    Returns
    ----------
    output : numpy array (float)
        the data in a 1D array    
    """
    size = os.path.getsize(fname)
    b = array.array('d')
    f = open(fname, "rb")
    b.fromfile(f, size // 8)
    f.close()
    lst = b.tolist()
    output = np.array(lst)

    if swapbyteorder == 1:
        output = output.newbyteorder()
    return output


def read_correlation_single(fname, swapbyteorder=0):
    """Read a single precision binary file without header information
      Does not reshape the data.
        
    Parameters
    ----------
    fname : string
        name of the double precision binary file
    
    swapbyteorder : integer
        1 (0) - do (not) swap byte order of each number
                may be required when reading/writing from different operating systems

    Returns
    ----------
    output : numpy array (float)
        the data in a 1D array    
    """
    size = os.path.getsize(fname)
    # print "sizes:", size, size/8
    b = array.array('f')
    f = open(fname, "r")
    b.fromfile(f, size // 4)
    f.close()
    lst = b.tolist()
    output = np.array(lst)

    if swapbyteorder == 1:
        output = output.newbyteorder()
    return output


def write_dbin(fname, data):
    """Write a double precision binary file without header information
        
    Parameters
    ----------
    fname : string
        name of the double precision binary file
    
    data : numpy array (float)
        data to be written to the file
    """    
    f = open(fname, "wb")
    fmt = '<' + 'd' * data.size
    bin_in = struct.pack(fmt, *data.flatten()[:])
    f.write(bin_in)
    f.close()


def h5read(filename, field="data/data1"):    
    """Read a image or array data from a h5 file 
        
    Parameters
    ----------
    filename : string
        name of the h5 file
    
    field : string
        field specifying the location of the data in the h5 file

    Returns
    ----------
    image : numpy array (float)
        image data    
    """
    h5file = h5py.File(filename, "r")
    # print field
    h5data = h5file[field]
    image = h5data[...]
    h5file.close()
    return image


def tifread(filename):
    """Read array data from a tif file
      Can read data greater than 2 dimensions 
        
    Parameters
    ----------
    filename : string
        name of the h5 file

    Returns
    ----------
    data : numpy array (float)
        image data    
    """
    data = np.asarray(Image.open(filename))
    # Checks if the array is 3D then flattens it
    if len(np.shape(data)) > 2:
        dig_data = np.ones((data.shape[0], data.shape[1]))
        xlen = len(data[:, 0, 0])
        ylen = len(data[0, :, 0])
        for x in range(0, xlen):
            for y in range(0, ylen):
                dig_data[x][y] = data[x][y][0]
        return dig_data
    else:
        return data


def read_image(fname, swapbyteorder=0, nx=-1, ny=-1, h5field="/data/data1",
                ser_index=0):
    """Read image array data based on file name extension
      Formats supported: dbin, bin, h5, tif, npy
        
    Parameters
    ----------
    fname : string
        name of the h5 file

    swapbyteorder=0 : integer 
        1 (0) - do (not) swap byte order of each number
        Only used for dbin or bin formats
        may be required when reading/writing from different operating systems
    
    nx (ny) : integer
        number of rows (columns) in the 2D array. 
        Only used for dbin and bin formats 

    h5field : string
        location of the image data in the h5 file

    ser_index : integer (for a currenlty unsupported file format)

    Returns
    ----------
    image : numpy array (float)
        image data    
    """
    fname = fname.rstrip('\n')
    fname = fname.rstrip('\r\n')
    print("DEBUG <io.py> fname", fname)
    if fname[-5:] == ".dbin":
        image = read_dbin(fname, swapbyteorder, nx=nx, ny=ny)
    elif fname[-4:] == ".bin":
        image = read_bin(fname, swapbyteorder).astype(np.float64)
#    elif fname[-4:] == ".ser":
#        image = read_ser(fname, ser_index)
    elif fname[-3:] == ".h5":
        image = h5read(fname, h5field)
    elif fname[-4:] == ".npy":
        image = np.load(fname)
    elif fname[-4:] == ".tif":  # Jack 28/10
        image = tifread(fname)  # Jack 28/10
    else:
        try:
            image = imageio.imread(fname)  # Jack 28/10 - this line is deprecated
        except:
            print("error reading image - unknown file location or extension. fname = ", fname)

    return image.astype(np.float64)

def write_image(fname, data, swapbyteorder=0, h5field="/data/data1"):
    """Write image array data based on file name extension
      Formats supported: dbin, npy
        
    Parameters
    ----------
    fname : string
        name of the h5 file

    swapbyteorder=0 : integer 
        1 (0) - do (not) swap byte order of each number
        Only used for dbin or bin formats
        may be required when reading/writing from different operating systems
    
    nx (ny) : integer
        number of rows (columns) in the 2D array. 
        Only used for dbin and bin formats 

    h5field : string
        location of the image data in the h5 file
    """
    fname = fname.rstrip('\n')
    fname = fname.rstrip('\r\n')
    if fname[-5:] == ".dbin":
        write_dbin(fname, data)
    elif fname[-4:] == ".npy":
        np.save(fname, data)
    else:
        print("<write_image()> : output extension not recognised - file not written")


def get_array_dimensions_from_file(fname):
    """
        
    Parameters
    ----------
    fname : string
        name of the image file 
        see read_image() for supported file types
        [note: currently doesn't support read_image() keywords]
    
    Returns
    ----------
    s : tuple
        shape of the image array

    """    
    image = read_image(fname)
    s = image.shape()
    return s

"""
def read_ser(fname, ser_index=0):
    data, calp, adl, sh = ser.read_ser(fname, [ser_index])
    image = data.dlist[0]
    return image


def read_ser_npatterns(fname):
    data, calp, adl, sh = ser.read_ser(fname, [1])
    return sh.ValNumElem


def read_ser_arraysize(fname):
    data, calp, adl, sh = ser.read_ser(fname, [1])
#    print "DEBUG <io.py; read_ser_arraysize> arraysize", calp.ArraysizeX, data.dlist[0].shape
#    return calp.ArraysizeX
"""

