# pypadf
***
Compute the pair-angle distribution function (PADF) from a fluctuation scattering dataset,
or
compute the pair-angle distribution function (PADF) from an atomic model (directly, without a diffraction calculation).


### Install

Clone this repo:

    git clone https://github.com/amartinrmit/pypadf.git

### Requirements

The following packages are required by pypadf. All testing has used the following versions, but later versions should also work.

    python==3.9
    numpy==1.26.3
    scipy==1.11.4
    matplotlib==3.8.0
    numba==0.59.0
    h5py==3.9.0
    imageio==2.33.1
    tqdm==4.67.1
    pyshtools==4.13.1


To install with conda (suggested), create and activate a conda environment and install the packages. This may take 10-15 minutes to install all from scratch.

    conda create --name pypadf python==3.9 -y
    conda activate pypadf
    conda install numpy==1.26.3 -y
    conda install scipy==1.11.4 -y
    conda install matplotlib==3.8.0 -y
    conda install numba==0.59.0 -y
    conda install h5py==3.9.0 -y
    conda install imageio==2.33.1 -y
    conda install tqdm==4.67.1 -y
    conda install conda-forge::pyshtools==4.13.1 -y

### Getting Started

For demonstration of the package, see instructions in the `demo` directory. This code has been run on Windows, Mac and Linux systems. We have not identified minimum OS requirements - if you can install the dependent python packages, it should work to our knowledge.

If you have a dataset, you can also implement the analysis pipeline via a jupyter notebook: padf_all_analysis.ipynb  
(By default this notebook can be run on the diffraction data generated in the demo script ./demo/hextest.py; run this demo first if you want to test the notebook.)
To implement the notebook on your dataset, you will need to copy the config files from demo and update the experimental parameters in the config files to match your data.
