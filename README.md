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

To run scripts individually, edit the hextest.py script and set the following parameters:
    
    test_diff = True
    test_corr = True
    test_mask = True
    test_padf = True
    test_plot = True
    test_dc   = True

Setting any of the following variables to `False` will skip that step in the pipeline, which are:

- `test_diff`: Diffraction simulation. Create and save a series of diffraction patterns of a hexagon structure.
- `test_corr`: Correlation caluclation. Calculate the scattering correlation function of the simulated diffraction patterns.
- `test_mask`: Mask calculation. Mask the correlation functions according to the masked regions of the detector.
- `test_padf`: PADF caclulation. Calculate the PADF from the scattering correlation function.
- `test_plot`: Plot various results.
- `test_dc`: Diffract and correlation. Simulates the diffraction pattern and calculates the correlation function, without saving the diffraction pattern inbetween (saving storage space)



