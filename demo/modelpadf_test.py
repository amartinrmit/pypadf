"""
A test script for the calculating the PADF from an atomic model
This calculation does not use the diffraction code
"""

import subprocess
import os
from pathlib import Path
def pstr( path, append_text=""):
    return str( (path/append_text).resolve())

configpath = Path( "./configs")
configfile = pstr( configpath, 'config_modelpadf.txt')

outpath = Path( "./output/modelpadf")
if not os.path.exists(outpath): os.makedirs(outpath)

sample =  "fcc_Si_3_3_3"
subjectxyz = pstr( configpath, f"{sample}.xyz")
method = "spharmonic"
mode = 'stm'
tag = f"{sample}_{method}_{mode}"

flags = f"--subjectxyz {subjectxyz} --outpath {pstr(outpath)} --tag {tag} --method {method} --mode {mode}  --tag {tag}"
result = subprocess.run(f'python ../modelpadf.py -c {configfile} '+flags,shell=True)


#
# plot the output
#
eo = "evens"  #or "odds"
padfname = pstr(outpath, tag+f"_mPADF_{eo}_sum.npy")
os.system(f"python ../plotfxs3d.py --config ./configs/config_sph_plot.txt --fname {padfname} --outpath {outpath} --power -1")

