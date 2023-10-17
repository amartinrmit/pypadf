# pypadf

Compute the pair-angle distribution function (PADF) from a fluctuation scattering dataset

## Build and Install

Clone this repo and move into it:

    git clone https://github.com/amartinrmit/pypadf.git
    cd pypadf

All of the following commands are assumed to be run from this directory.

Create directory to save outputs. Template config files will save data to these directories, but can be edited to a directory of your choosing.

    mkdir ./tmp
    mkdir ./tmp/diff
    mkdir ./tmp/mask
    mkdir ./tmp/corr

Requirements can be installed with conda or pip. To install with conda (suggested), create and activate a conda environment and install the packages.

    conda create --name pypadf python==3.8
    conda activate pypadf
    conda install numpy -y
    conda install scipy -y
    conda install matplotlib -y
    conda install numba -y
    conda install h5py -y
    conda install imageio -y

Alternatively, the same packages can be installed via pip (untested).

    pip install numpy
    pip install scipy
    pip install matplotlib
    pip install numba
    pip install h5py
    pip install imageio
    

## Worked example

#### Step 1: Simulate Diffraction Patterns

First step will be to simulate a some diffraction patterns. Parameters can be read from a config file.

    python diffract.py --config ./configs/config_hex_diff.txt

This will create 2 diffraction patterns and save them to `./tmp/dp/`. 
Any parameter in the config file can be overided on the commandline. To generate more patterns:

    python diffract.py --config ./configs/config_hex_diff.txt --npatterns 1000

To see all options:
    
    python diffract.py --help




#### Step 2: Inspect Diffraction Pattern

To inspect a diffraction pattern:

    python plotdiffraction.py --fname ./tmp/diff/hex_1.npy





#### Step 2.5: Make mask

We will create a mask file that is 1 for every pixel in the difraction pattern (essentially no mask) and save it to `./tmp/mask`.

    python3 make-mask.py

#### Step 3: Correlate 

Correlate 2 diffraction patterns.

    python difftocorr.py --config ./configs/config_hex_corr.txt

Again, parameters can be overridden on the command line.

    python difftocorr.py --config ./configs/config_hex_corr.txt --npatterns 1000


This will generate new config files to plot the correlation and create the PADF

#### View correlation

The generated config file from running `difftocorr.py` can be used to plot the q1=q2 plane of the correlation function.

    python plotfxs3d.py --config ./configs/config_hex_a_plot.txt

To better see the correlation intensity, try chaning the colorscale:

    python plotfxs3d.py --config ./configs/config_hex_a_plot.txt --chigh 0.00005 --clow -0.00005

You can also try correlating with fewer patterns, and replotting to see the difference. 


#### PADF

The generated config file from running `difftocorr.py` can be used to generate the padf. This will save the output function to the same fold as the input correlation file.

    python corrtopadf.py --config ./configs/config_hex_a_padf.txt



# Notes:
- is there a way to set up the configs so we dont require a mask?
- is there a way to turn this off auto-generation of configs and just include template configs for the padf and plotting steps?
- if we want to keep the auto-generation, we should be more selective in the .gitignore to remove the pregenerated ones?
- Plotting for padf? when I tried running the script with the same config but chaning the --fname, it just plotted the q correlation again.



