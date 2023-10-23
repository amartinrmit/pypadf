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

<<<<<<< HEAD
### All-In-One Demonstration

Move to the `demo` directory, and run the `hextest.py` script to run through each step of the PADF calculation.

    cd demo
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




### Individual Script Demonstration (Linux)

The following steps go into further depth of running each step of the PADF calculation.


#### Create Output Directories 
To illustrate the pypadf package, we will run the scripts with provided config files in `./demo/configs`. These template config files save certain outputs to directories which we will now create. From the top level directory, move into the demo folder and create the output directories.

    cd demo
    mkdir ./output
    mkdir ./output/diff
    mkdir ./output/mask
    mkdir ./output/corr
    mkdir ./output/dc

Alternatively, the config files can be edited for an output directory of your choosing.

#### Simulate Diffraction Patterns

We will now simulate some diffraction patterns. Parameters are read from a config file.

    python ../diffract.py --config ./configs/config_hex_diff.txt

This will create 6 diffraction patterns and save them to the output directory `./demo/output/diff`. 

To see all options:
    
    python ../diffract.py --help


#### Inspect Diffraction Pattern

To inspect a diffraction pattern:

    python ../plotdiffraction.py --fname ./output/diff/hex_0.npy

##### Make Mask

We will create a mask file that is 1 for every pixel in the difraction pattern (essentially no mask) and save it to `./demo/output/mask`.

    python ../make-mask.py ./output/diff/hex_0.npy ./output/mask/hex_mask.npy

#### Correlate the Diffraction Patterns

Correlate 6 diffraction patterns. The number of patterns will be split into two correlation functions, an A half from 3 patterns, and a B half from the other 3 patterns.

    python ../difftocorr.py --config ./configs/config_hex_corr.txt


This will generate new config files to plot the correlation and create the PADF

#### Diffraction and correlate

If you don't want to save 1000 diffraction patterns, then you can run diffract and correlate to generate the patterns, and correlate directly.

    python ../diffract_and_correlate.py --config ./configs/config_hex_dc.txt

#### View correlation

The generated config file from running `difftocorr.py` can be used to plot the q1=q2 plane of the correlation function.

    python ../plotfxs3d.py --config ./configs/config_hex_a_corr_plot.txt

To better see the correlation intensity, try chaning the colorscale:

    python ../plotfxs3d.py --config ./configs/config_hex_a_corr_plot.txt --chigh 0.00005 --clow -0.00005

You can also try correlating with fewer patterns, and replotting to see the difference. 


#### PADF

The generated config file from running `difftocorr.py` can be used to generate the PADF. This will save the output function to the same fold as the input correlation file.

    python ../corrtopadf.py --config ./configs/config_hex_a_padf.txt
=======

### Getting Started
>>>>>>> 4a66775 (split readme in demo)

For demonstration of the package, see instructions in the `demo` directory.
