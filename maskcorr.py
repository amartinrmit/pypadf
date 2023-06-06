
import numpy as np
import matplotlib.pyplot as plt
import params.paramsMASK as params
import os

#
# set up the parameter class
#
p = params.paramsMASK()

print("maskcorr.py")
print("pypadf version ",p.version)
#
# read input parameters from file / command line
#
p.read_config_file()

#
# load the correlation file
#
corrfile = p.path_to_string( p.corrfile )
corr = np.load( corrfile )
print( "Correlation loaded:", corrfile ) 

s = corr.shape


#
# mask theta
#
if p.maskth:
     ithw = int(s[2]*p.thlim/360)
     corr[:,:,:ithw]  = 0
     corr[:,:,-ithw:] = 0


#print("27 DEBUG min max", np.min(corr), np.max(corr))
#
# subtract the mean value
#
if p.submean:
        if p.thlimnorm>0:
            ith0 = int(s[2]*p.thlimnorm/360)
        else:
            ith0 = 0

        for i in range(s[0]):
            for j in range(s[1]):
                corr[i,j,:] += -np.average(corr[i,j,ith0:-ith0])

#
# sintheta correction
#
if p.sintheta:
        nth = corr.shape[2]
        th = 2*np.pi*np.arange(nth)/nth
        th += th[1]/2
        sth = np.abs(np.sin(th))
        for i in range(s[0]):
            for j in range(s[1]):
                corr[i,j,:] *= sth

#
# mask q
#
if p.maskq:
        iqhigh = int( s[0]*(p.qmaskhigh-p.qmin)/(p.qmax-p.qmin))
        if iqhigh<0 or iqhigh>=s[0]:
              print("qmaskhigh value is out of (qmin,qmax) range")
              exit()

        iqlow = int( s[0]*(p.qmasklow-p.qmin)/(p.qmax-p.qmin))
        if iqlow<0 or iqlow>=s[0]:
              print("qmasklow value is out of (qmin,qmax) range")
              exit()

        corr[:iqlow,:iqlow,:] = 0
        corr[iqhigh:,iqhigh:,:] = 0

#
# mask theta
#
if p.maskth:
     ithw = int(s[2]*p.thlim/360)
     corr[:,:,:ithw]  = 0
     corr[:,:,-ithw:] = 0

#
# check for nans
#
corr[np.isnan(corr)] = 0.0

#print("DEBUG min max", np.min(corr), np.max(corr), iqlow, iqhigh, s[0])

#
# output
#
#outpath = p.path_to_string(outpath)
outname = p.makefname( p.outpath, os.path.basename(p.corrfile)[:-4],p.suffix,".npy") 
print( "Output file: ", outname ) 
np.save( outname, corr)
