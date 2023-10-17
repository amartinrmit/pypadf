# pypadf

Compute the pair-angle distribution function (PADF) from a fluctuation scattering dataset

## Build and Install

Clone this repo and move into it:

    git clone https://github.com/amartinrmit/pypadf.git
    cd pypadf

All of the following commands are assumed to be run from this directory.

Create directory to save outputs. Template config files will save data to these directories, but can be edited to a directory of your choosing.

    mkdir ./demo/output
    mkdir ./demo/output/diff
    mkdir ./demo/output/mask
    mkdir ./demo/output/corr
    mkdir ./demo/output/padf

Requirements can be installed with conda or pip. To install with conda (suggested), create and activate a conda environment and install the packages.

    conda create --name pypadf python==3.8
    conda activate pypadf
    conda install numpy -y
    conda install scipy -y
    conda install matplotlib -y
    conda install numba -y
    conda install h5py -y
    conda install imageio -y

## Worked example

#### Step 1: Simulate Diffraction Patterns

First step will be to simulate a some diffraction patterns. Parameters can be read from a config file.

    python diffract.py --config ./demo/configs/config_hex_diff.txt

This will create 2 diffraction patterns and save them to `./demo/output/diff`. 
Any parameter in the config file can be overided on the commandline. To generate more patterns:

    python diffract.py --config ./demo/configs/config_hex_diff.txt --npatterns 1000

To see all options:
    
    python diffract.py --help




#### Step 2: Inspect Diffraction Pattern

To inspect a diffraction pattern:

    python plotdiffraction.py --fname ./demo/output/diff/hex_0.npy





#### Step 2.5: Make mask

We will create a mask file that is 1 for every pixel in the difraction pattern (essentially no mask) and save it to `./demo/output/mask`.

    python make-mask.py ./demo/output/diff/hex_0.npy ./demo/output/mask/hex_mask.npy

#### Step 3: Correlate 

Correlate 2 diffraction patterns.

    python difftocorr.py --config ./demo/configs/config_hex_corr.txt

Again, parameters can be overridden on the command line.

    python difftocorr.py --config ./demo/configs/config_hex_corr.txt --npatterns 1000


This will generate new config files to plot the correlation and create the PADF

#### View correlation

The generated config file from running `difftocorr.py` can be used to plot the q1=q2 plane of the correlation function.

    python plotfxs3d.py --config ./demo/configs/config_hex_a_plot.txt

To better see the correlation intensity, try chaning the colorscale:

    python plotfxs3d.py --config ./demo/configs/config_hex_a_plot.txt --chigh 0.00005 --clow -0.00005

You can also try correlating with fewer patterns, and replotting to see the difference. 


#### PADF

The generated config file from running `difftocorr.py` can be used to generate the padf. This will save the output function to the same fold as the input correlation file.

    python corrtopadf.py --config ./demo/configs/config_hex_a_padf.txt



# Notes:
- I dont get why a new file hexagon_shifted.pdb gets created, it doesn't fix the runtime warning for invalid value encountered in scalar remainder during diffraction pattern generation

- when generating 1000 diffraction patterns, it says x/1000 patterns made
- when generating 1000 correaltion functions, it says x/999 patterns made

- plotfxs3d: there is a print statement that says "section extracted: reqr". this probably shouldn't be there?

- corrtopadf: step one says correlation to blrr but should be blqq








- is there a way to set up the configs so we dont require a mask?
- is there a way to turn this off auto-generation of configs and just include template configs for the padf and plotting steps?
- if we want to keep the auto-generation, we should be more selective in the .gitignore to remove the pregenerated ones?
- Plotting for padf? when I tried running the script with the same config but chaning the --fname, it just plotted the q correlation again.



