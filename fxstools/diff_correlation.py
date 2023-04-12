
import numpy as np
import matplotlib.pyplot as plt

rmax = 3.27*2

path = "/Users/andrewmartin/cloudstor/Work/Research/SF/codes/code_projects/py3padf_0.3/dev/"
fname = "corr_predicted.npy"
corr = np.load( path+fname )
s = corr.shape
path2 = path
fname2 = "twocarbons_test1_correlation_sum.npy"
corr2 = np.load( path2+fname2 )

norm1 = np.sqrt( np.sum( corr**2))
norm2 = np.sqrt( np.sum( corr2**2))

corr2 *= np.abs(np.max(corr[25:,25:,:])/np.max(corr2[25:,25:,:]))   #norm1/norm2
difference = np.sqrt( np.sum( (corr-corr2)**2))
print( "Correlation difference :", difference, norm1, norm2 )

"""
for i in np.arange(corr.shape[0]):
    for j in np.arange(corr.shape[0]):
        if (i !=0)and(j!=0):
            corr[i,j,:] *= (j*i)**2 
            corr2[i,j,:] *= (j*i)**2
"""

rvals = rmax*np.arange(corr.shape[0])/corr.shape[0]

disp = np.zeros( (corr.shape[0], corr.shape[2] ) )
disp2 = np.zeros( (corr2.shape[0], corr2.shape[2] ) )
disp3 = np.zeros( (corr2.shape[0], corr2.shape[2] ) )
for i in np.arange(corr.shape[0]):
    disp[i,:] = corr[i,i,:]
    disp2[i,:] = corr2[i,i,:]
    disp3[i,:] = corr[i,i,:] - corr2[i,i,:]
#disp = corr[:,:,0]

xt = [0,90,180,270]
scl, sc = 1.0, 1.0
x0,x1 = 0, s[0]
print(x0,x1)
r0 = x0*rmax/s[0]
plt.xticks(xt)
r1 = x1*rmax/s[0]
plt.subplot(131)
plt.imshow(disp[x0:x1,:], extent=[0,360,r0,r1], origin='lower',aspect=100)
plt.clim([np.min(disp)*scl,np.max(disp)*sc])
plt.colorbar(fraction=0.046, pad=0.04)
plt.title("predicted from the atom padf")
plt.xticks(xt)

plt.subplot(132)
plt.imshow(disp2[x0:x1,:],  extent=[0,360,r0,r1], origin='lower',aspect=100)
plt.clim([np.min(disp2)*scl,np.max(disp2)*sc])
plt.colorbar(fraction=0.046, pad=0.04)
plt.title("from diffraction code")
plt.xticks(xt)

plt.subplot(133)
plt.imshow(disp3[x0:x1,:],  extent=[0,360,r0,r1], origin='lower',aspect=100)
plt.clim([np.min(disp)*scl,np.max(disp)*sc])
plt.colorbar(fraction=0.046, pad=0.04)
plt.title("difference")
plt.xticks(xt)
plt.tight_layout()

plt.figure()
line = 8
fact = np.max(disp[line,:])/np.max(disp2[line,:])
plt.plot( disp[line,:] )
plt.plot( disp2[line,:]*fact )
plt.draw()
plt.show()
