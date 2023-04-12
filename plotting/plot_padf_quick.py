
import numpy as np
import matplotlib.pyplot as plt

data = np.load("/Users/E38496/Work/Research/SF/codes/code_projects/py3padf_0.4/test/corr/testdp_padf.npy")
#data = np.load("twopoints_test1_correlation_sum.npy")
#data = np.load("corr_interp2_DEBUG.npy")

print(data.shape, np.min(data), np.max(data) )

disp = np.zeros( (data.shape[0],data.shape[2]))


for i in np.arange(data.shape[0]):
    disp[i,:] = data[i,i,:] - np.average(data[i,i,:])
    if i>0:
      disp[i,:] *= (i*i)**2


rmax = 6.0 # 65.34

r0,r1 = 0,data.shape[0]-0
rmindisp = r0*rmax/data.shape[0]
rmaxdisp = r1*rmax/data.shape[0]

scl, sc = 0, 1 
plt.imshow( disp[r0:r1,:], origin='lower', extent=[0,360,rmindisp,rmaxdisp], aspect=25) 
plt.clim([-np.max(disp[r0:r1,:])*scl,np.max(disp[r0:r1,:])*sc])

plt.figure()
plt.plot( disp[r0:r1,data.shape[2]//2])

rdata = np.arange(r1-r0)*(rmaxdisp-rmindisp)/(r1-r0) + rmindisp
#np.save( "th180.npy", np.array([rdata,disp[r0:r1,data.shape[2]//2]]))

plt.show()
