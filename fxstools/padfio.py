"""
padfio.py

input and output functions for correlation and padf calculations

Functions:
    read_dbin
    read_bin
    read_correlation
    read_correlation_single
    write_dbin
    h5read
    tifread
    read_image
    get_array_dimensions_from_file
"""

import array
import imageio
import os
import struct
import h5py
from PIL import Image
import pathlib

#import fxstools.amser as ser
import numpy as np


def read_dbin(fname, swapbyteorder=0, nx=-1, ny=-1):
    """
    Read a double precision binary file

    Parameters
    ----------
    fname : str
        name of the input image file (including path)
        supported file formats: .dbin, .bin, .npy,  .h5, .tif
    
    nx : int
        number of rows in input image array
        (required for binary file reading .bin, .dbin)

    ny : int
        number of columns in input image array
        (required for binary file reading .bin, .dbin)
    
    swapbyteorder : int
        for binary file reading across file systems
        0 - don't swap byte order in bin, dbin file
        1 - swapbyte order in bin, dbin file

    Returns
    -------
    output : numpy array (float)
        2D array (floats) containing the data

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
    """
    Read a single precision binary file

    Parameters
    ----------
    fname : str
        name of the input image file (including path)
        supported file formats: .dbin, .bin, .npy,  .h5, .tif
    
    nx : int
        number of rows in input image array
        (required for binary file reading .bin, .dbin)

    ny : int
        number of columns in input image array
        (required for binary file reading .bin, .dbin)
    
    swapbyteorder : int
        for binary file reading across file systems
        0 - don't swap byte order in bin, dbin file
        1 - swapbyte order in bin, dbin file

    Returns
    -------
    output : numpy array (float)
        2D array (floats) containing the data

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
    """
    Read a double precision binary file 
    without reshaping into 2D array

    Parameters
    ----------
    fname : str
        name of the input image file (including path)
        supported file formats: .dbin, .bin, .npy,  .h5, .tif
    
    swapbyteorder : int
        for binary file reading across file systems
        0 - don't swap byte order in bin, dbin file
        1 - swapbyte order in bin, dbin file

    Returns
    -------
    output : numpy array (float)
        1D array (floats) containing the data

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
    """
    Read a single precision binary file 
    without reshaping into 2D array

    Parameters
    ----------
    fname : str
        name of the input image file (including path)
        supported file formats: .dbin, .bin, .npy,  .h5, .tif
    
    swapbyteorder : int
        for binary file reading across file systems
        0 - don't swap byte order in bin, dbin file
        1 - swapbyte order in bin, dbin file

    Returns
    -------
    output : numpy array (float)
        1D array (floats) containing the data

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
    """
    Write a numpy array as double precision binary file

    Parameters
    ----------
    fname : str
        name of the output file (including path)

    Returns
    -------

    """
    f = open(fname, "wb")
    fmt = '<' + 'd' * data.size
    bin_in = struct.pack(fmt, *data.flatten()[:])
    f.write(bin_in)
    f.close()


def h5read(filename, field="data/data1"):
    """
    Read image (intensity) data from a .h5 file

    Parameters
    ----------
    filename : str
        name of the input h5 file (including path)

    Returns
    -------
    image : numpy array (float)
        array containing the image data
        intensity values only; no colour dimension

    """
    h5file = h5py.File(filename, "r")
    # print field
    h5data = h5file[field]
    image = h5data[...]
    h5file.close()
    return image


def tifread(filename):
    """
    Read a .tiff file containing floating point intensity values

    Parameters
    ----------
    filename : str
        name of the input tiff file (including path)

    Returns
    -------
    data : numpy array (float)
        image array with intensity values

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


def read_image(fname, swapbyteorder=0, nx=-1, ny=-1, h5field="/data/data1"):
    """
    Read an image file based on file extension

    Parameters
    ----------
    fname : str
        name of the input image file (including path)
        supported file formats: .dbin, .bin, .npy,  .h5, .tif
    
    nx : int
        number of rows in input image array
        (required for binary file reading .bin, .dbin)

    ny : int
        number of columns in input image array
        (required for binary file reading .bin, .dbin)
    
    swapbyteorder : int
        for binary file reading across file systems
        0 - don't swap byte order in bin, dbin file
        1 - swapbute order in bin, dbin file 

    h5field : str
        h5 field containing the array of image data
        only used for h5 file reading

    Returns
    -------
    image : numpy array (float)
        array containing the image data
        intensity values only; no colour dimension
    """
    fname = fname.rstrip('\n')
    fname = fname.rstrip('\r\n')
    #print("DEBUG <io.py> fname", fname)
    if fname[-5:] == ".dbin":
        image = read_dbin(fname, swapbyteorder, nx=nx, ny=ny)
    elif fname[-4:] == ".bin":
        image = read_bin(fname, swapbyteorder).astype(np.float64)
    elif fname[-4:] == ".npy":
        image = np.load(fname)
    #elif fname[-4:] == ".ser":
    #    image = read_ser(fname, ser_index)
    elif fname[-3:] == ".h5":
        image = h5read(fname, h5field)
    elif fname[-4:] == ".tif":  # Jack 28/10
        image = tifread(fname)  # Jack 28/10
    else:
        try:
            image = imageio.imread(fname)  # Jack 28/10 - this line is deprecated
        except:
            print("error reading image - unknown file location or extension. fname = ", fname)

    return image.astype(np.float64)


def makefname( path, tag, suffix, fext):
        outname = path / (tag+suffix+fext)
        outname = str(outname.resolve())
        return outname



def get_array_dimensions_from_file(fname):
    """
    Open one image file and return the array shape

    Parameters
    ----------
    fname : str
        image file name to read

    Returns
    -------
    s : tuple (int)
        shape of image array
    """
    image = read_image(fname)
    s = image.shape()
    return s
