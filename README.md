# pypadf
***
Compute the pair-angle distribution function (PADF) from a fluctuation scattering dataset.

### Install

Clone this repo:

    git clone https://github.com/amartinrmit/pypadf.git

### Requirements

The following packages are required by pypadf:

    numpy, scipy, matplotlib, numba, h5py, imageio

To install with conda (suggested), create and activate a conda environment and install the packages.

    conda create --name pypadf python==3.9 -y
    conda activate pypadf
    conda install numpy -y
    conda install scipy -y
    conda install matplotlib -y
    conda install numba -y
    conda install h5py -y
    conda install imageio -y

### Demonstration

Move to the `demo` directory, and run the `hextest.py` script to run through each step of the PADF calculation.

    python hextest.py

This will automattically make an output directory to save various simulated quantaties during the PADF calculation. If the `hextest.py` script runs with no errors, then everything should be installed correctly.


