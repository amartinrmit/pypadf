import numpy as np
import matplotlib.pyplot as plt


data = np.load("./fsing.npy")

s =  data.shape

# some test things
z = np.arange(100)*2/100 - 1
th = np.arccos( z )
ith = 100*th/(2*np.pi)

print( ith )

# actual plotting
for i in np.arange(s[0]):

    plt.plot( data[i,:] )

plt.draw()
plt.show()

mat = np.load( "./testmat.npy" )
normmat = np.dot( mat, np.transpose(mat) )

s = mat.shape
diag = np.zeros(s[0])
for i in np.arange(s[0]):
    i2 = 2*i
    diag[i] = (normmat[i,i]*(2*i2+1) - 0)  #/ (i*i*i) +1

for i in np.arange(s[0]):
    for j in np.arange(s[0]):
        i2, j2 = i*2, j*2
        normmat[i,j] *= np.sqrt(2*i2+1)*np.sqrt(2*j2+1)

plt.imshow( normmat )
plt.figure()
plt.plot( diag )
plt.ylim([0,np.max(diag)])
plt.draw()
plt.show()
