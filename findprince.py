'''find principle axes of a halo and their satellites from the distribution of satellite, following the fortran code written by Noam. I. Libeskind
'''
import numpy as np
from numpy import dot, sqrt
from helpers import coordwrap2

def findprince(x,y,z):
	#returns principle axes in terms of a 3x3 matrix, where each column is a principle axis in one direction
	
	#(x,y,z) original positions
	N = len(x)
	xx=0.
	yy=0.
	zz=0.

	xz=0.
	yz=0.
	xy=0.
	for i in range(N):
		xx += (x[i]**2e0)
		yy += (y[i]**2e0)
		zz += (z[i]**2e0)
     
		xy += (x[i]*y[i])
		yz += (y[i]*z[i])
		xz += (x[i]*z[i])
		
	Inert = np.zeros((3,3))
	
	#entering inertial tensor computed from the position of the particles, or satellite galaxies, in this case
	Inert[0][0] = xx
	Inert[0][1] = xy
	Inert[0][2] = xz
	Inert[1][0] = xy
	Inert[1][1] = yy
	Inert[1][2] = yz
	Inert[2][0] = xz
	Inert[2][1] = yz
	Inert[2][2] = zz
	
	#eigenvalues array w and eigenvectors matrix v
	w, v = np.linalg.eig(Inert)

	idx = w.argsort()[::-1]  #sort eigenvalues from large to small
	w = w[idx]
#	print "inertial tensor eigenvalues", w
	v = v[:,idx]	#eigenvectors correspond to eigenvalues from large to small
#	print v
	return v
	
def rotate(x, y, z, v):	
	#(nx,ny,nz) and rotated positions
#	for i in range(3):
#		v[i] =v[i]/np.linalg.norm(v[i])
	
	N = len(x)
	nx = np.zeros(N)
	ny = np.zeros(N)
	nz = np.zeros(N)
	
	for i in range(N):
#		[nx[i], ny[i], nz[i]] = dot(np.array([x[i],y[i],z[i]]), v).tolist()
		nx[i] = v[0][0]*x[i] + v[1][0]*y[i] + v[2][0]*z[i]
		ny[i] = v[0][1]*x[i] + v[1][1]*y[i] + v[2][1]*z[i]
		nz[i] = v[0][2]*x[i] + v[1][2]*y[i] + v[2][2]*z[i]
	
	return nx, ny, nz
