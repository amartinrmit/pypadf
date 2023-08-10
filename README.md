# pypadf

Compute the pair-angle distribution function (PADF) from a fluctuation scattering dataset





## Build

Clone this repo and move into it:

    git clone https://github.com/amartinrmit/pypadf.git
    cd pypadf

All of the following commands are assumed to be run from this directory.

Create directory to save outputs. Template config files will save data to these directories, but can be edited to a directory of your choosing.

    mkdir ./tmp
    mkdir ./tmp/diff
    mkdir ./tmp/mask
    mkdir ./tmp/corr
    mkdir ./tmp/padf


DEV: For a clean install, remove any previous pypadf environments:

    conda remove -n pypadf --all

Create and activate a conda environment. Not nessercary, but advised.

    conda create --name pypadf
    conda activate pypadf

Install required packages. The following are conda instructions, but pip should work as well.

    conda install numpy -y
    conda install scipy -y
    conda install matplotlib -y
    conda install numba -y
    conda install h5py -y
    conda install imageio -y



## Worked example



###  Simulate diffraction patterns:

First step will be to simulate a some diffraction patterns. Parameters can be read from a config file.

    python diffract.py --config ./configs/config_hex_dp.txt

This will create 2 diffraction patterns and save them to `./data/dp/`. 
Any parameter in the config file can be overided on the commandline. To generate more patterns:

    python diffract.py --config ./configs/config_hex_dp.txt --npatterns 1000

To see all options:
    
    python diffract.py --help




## view diffraction

To inspect a diffraction pattern: AM can this also take a config file? for consistancy? low priority

    python plot-dp.py -f ./data/dp/hex_1.npy

## Make mask
A mask is required. just takes the first hex pattern and makes it an array of 1s
AM I couldn't work out how to make the configs run without a mask

    python3 make-mask.py

## correlate

Correlate 2 diffraction patterns.

    python dp-to-corr.py --config ./configs/config_hex_corr.txt

Again, parameters can be overridden on the command line.

    python dp-to-corr.py --config ./configs/config_hex_corr.txt --npatterns 10000

## view correlation

    python plotfxs3d.py --config ./configs/config_hex_a_plot.txt
    python plotfxs3d.py --config ./configs/config_hex_a_plot.txt --chigh 0.00005 --clow -0.00005


## 



















