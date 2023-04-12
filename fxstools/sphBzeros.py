### recursive method: computes zeros ranges of Jn(r,n) from zeros of Jn(r,n-1)
### (also for zeros of (rJn(r,n))')
### pros : you are certain to find the right zeros values;
### cons : all zeros of the n-1 previous Jn have to be computed;
### note : Jn(r,0) = sin(r)/r

import numpy as np
from scipy import arange, pi, sqrt, zeros
from scipy.special import jv, jvp, yv
from scipy.optimize import brentq
from sys import argv
from pylab import *

def Jn(r,n):
  return (sqrt(pi/(2*r))*jv(n+0.5,r))

def Jn_zeros(n,nt):
  zerosj = zeros((n+1, nt), dtype=np.float)
  zerosj[0] = arange(1,nt+1)*pi
  points = arange(1,nt+n+1)*pi
  racines = zeros(nt+n, dtype=np.float)
  for i in range(1,n+1):
    for j in range(nt+n-i):
      foo, r = brentq(Jn, points[j], points[j+1], (i,), full_output=True)
      racines[j] = foo
      #if j == 0: print( i, j, points[j], points[j+1], foo, r.converged)
    points = racines
    zerosj[i][:nt] = racines[:nt]
  return (zerosj)

def rJnp(r,n):
  return (0.5*sqrt(pi/(2*r))*jv(n+0.5,r) + sqrt(pi*r/2)*jvp(n+0.5,r))

def rJnp_zeros(n,nt):
  zerosj = zeros((n+1, nt), dtype=np.float)
  zerosj[0] = (2.*arange(1,nt+1)-1)*pi/2
  points = (2.*arange(1,nt+n+1)-1)*pi/2
  racines = zeros(nt+n, dtype=np.float)
  for i in range(1,n+1):
    for j in range(nt+n-i):
      foo = brentq(rJnp, points[j], points[j+1], (i,))
      racines[j] = foo
    points = racines
    zerosj[i][:nt] = racines[:nt]
  return (zerosj)
