# All-In-One Demonstration

Run the `hextest.py` script to run through all steps of the PADF calculation.

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




### Worked Example (Linux)

Alternatively, you can run each script individually for greater control of the analysis pipeline. The following steps go into further depth of running each step of the PADF calculation individually.


#### Create Output Directories 
To illustrate the pypadf package, we will run the scripts with provided config files in `./demo/configs`. These template config files save certain outputs to directories which we will now create.

    mkdir -p ./output
    mkdir -p ./output/diff
    mkdir -p ./output/mask
    mkdir -p ./output/corr
    mkdir -p ./output/dc
    mkdir -p ./output/padf


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

#### Correlate the Diffraction Patterns

Correlate 6 diffraction patterns. The number of patterns will be split into two correlation functions, an A half from 3 patterns (`output/corr/hex_a_correlation_sum.npy`) and a B half from the other 3 patterns (`output/corr/hex_a_correlation_sum.npy`).

    python ../difftocorr.py --config ./configs/config_hex_corr.txt


#### Diffraction and correlate

If you don't want to save 1000 diffraction patterns, then you can run diffract and correlate to generate the patterns, and correlate directly.

    python ../diffract_and_correlate.py --config ./configs/config_hex_dc.txt


#### PADF

The generated config file from running `difftocorr.py` can be used to generate the PADF. This will save the output function to the same fold as the input correlation file.

    python ../corrtopadf.py --config ./configs/config_hex_padf.txt



#### View correlation (Not working)

The generated config file from running `difftocorr.py` can be used to plot the q1=q2 plane of the correlation function.

    python ../plotfxs3d.py --config ./configs/config_hex_plot.txt







##### Make Mask (Not needed)

We will create a mask file that is 1 for every pixel in the difraction pattern (essentially no mask) and save it to `./demo/output/mask`.

    python ../make-mask.py ./output/diff/hex_0.npy ./output/mask/hex_mask.npy



