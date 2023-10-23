# All-In-One Demonstration

Run the `hextest.py` script to run through all steps of the PADF calculation.

    python hextest.py

This will automattically make an output directory to save various simulated quantaties during the PADF calculation. If the `hextest.py` script runs with no errors, then everything should be installed correctly.


To run scripts individually, edit the `hextest.py` script and set the following parameters:
    
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




# Worked Example (Linux)

Alternatively, you can run each script individually for greater control of the analysis pipeline. The following steps go into further depth of running each step of the PADF calculation individually.


## 1) Create Output Directories 
To illustrate the pypadf package, we will run the scripts with provided config files in `./demo/configs`. These template config files save certain outputs to directories which we will now create.

    mkdir -p ./output
    mkdir -p ./output/diff
    mkdir -p ./output/mask
    mkdir -p ./output/corr
    mkdir -p ./output/dc
    mkdir -p ./output/padf
    mkdir -p ./output/figs


Alternatively, the config files can be edited for an output directory of your choosing.

## 2) Simulate Diffraction Patterns

We will now simulate some diffraction patterns with the `diffract.py` script. Parameters are read from a config file.

    python ../diffract.py --config ./configs/config_hex_diff.txt

This will create 6 diffraction patterns and save them to the output directory `./demo/output/diff`. 

To see all options:
    
    python ../diffract.py --help


## 3) Inspect A Diffraction Pattern

To inspect a diffraction pattern, use the script `plotdiffraction.py`:

    python ../plotdiffraction.py --fname ./output/diff/hex_0.npy

## 4) Correlate the Diffraction Patterns

Once we have generated a set of diffraction patterns, we can correlate them with the `diftocorr.py` script. Of the 6 patterns we generated in the previous step, they will be split into two correlation functions. An A half from 3 patterns (`./output/corr/hex_a_correlation_sum.npy`) and a B half from the other 3 patterns (`./output/corr/hex_b_correlation_sum.npy`).

    python ../difftocorr.py --config ./configs/config_hex_corr.txt


Alternatively, we have supplied the script `diffract_and_correlate.py` that will simulate a set of diffraction patterns and correlate them directly. This will not save the itermeadite diffractiopn patterns, and can save storage space. 

    python ../diffract_and_correlate.py --config ./configs/config_hex_dc.txt


## 5) PADF

The script `corrtopadf.py` will generate the PADF from the correlation function made in the previous step, specifically the A half.

    python ../corrtopadf.py --config ./configs/config_hex_padf.txt


## 6) View PADF

The script `plotfx23d.py` will run produce a plot of r1=r2 plane within the PADF we generated in the previous steps.

    python ../plotfxs3d.py --config ./configs/config_hex_plot.txt





