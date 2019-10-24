import numpy as np
from numpy import random
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from findprince import *
from helpers import *
import sys
from numpy.linalg import inv
from constant import *

def plotlocal(sathalos, plotname):
	plt.clf()
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	[satofAx, satofAy, satofAz, satofBx, satofBy, satofBz] = sathalos
	
	for i in range(1000):
		ax.scatter(satofAx[i], satofAy[i], satofAz[i], zdir='z', color = 'black')

	for i in range(1000):
		ax.scatter(satofBx[i], satofBy[i], satofBz[i], zdir='z', color = 'blue')
	
	ax.set_xlim(-10,10)
	ax.set_ylim(-10,10)
	ax.set_zlim(-10,10)	
	fig.savefig(plotname, dpi=500)
	plt.close()
	

if __name__ == "__main__":
	A = np.array([[63.5, 20, 5], [0.4, 21, 4.6],[43,5,6], [63.5, 20, 5], [0.4, 21, 4.6],[43,5,6]])
	B = np.array([[1,0,1],[0,0,1],[0,1,0], [63.5, 20, 5], [0.4, 21, 4.6],[43,5,6]])
	b = 63.5
	a = 1.5
#	A = np.array([ 29.66408,  45.73954,  32.61423])
#	B = np.array([ 29.54962,  45.36013,  32.92216])
#	C = np.array([ 30.00162,  45.50572,  32.61423])
#	print(coorddistance(A, B))
#	print(coorddistance(middlecoordmpc(A,B), C))
#	boxmin = coordwrap2(A - 1.5e0 * LGdmax_sqrt) #within a box extended by 1.5Mpc * 1.5 in three direction 
#	boxmax = coordwrap2(A + 1.5e0 * LGdmax_sqrt) #with the jth halo at its center	
#	print(np.ndarray.all(np.mod(C-boxmin-dim/2e0, dim)-dim/2e0 > 0 &  np.ndarray.all(np.mod(C-boxmax-dim/2e0, dim)-dim/2e0 < 0)))
#	print(np.mod(C-boxmin-dim/2e0, dim)-dim/2e0, np.mod(C-boxmax-dim/2e0, dim)-dim/2e0, np.ndarray.any(C!=A))
	
	coord1 = np.array([np.array([63.5, 33.3, 59]), np.array([63.5, 20, 5]), np.array([0.4, 1, 4.6])])
	coord2 = np.array([np.array([0., 63.7, 58]), np.array([0.4, 21, 4.6]), np.array([4,34,6])])
	
#	print(coord1 - np.array([1,2,1]))
	res = np.linalg.norm(B, axis =1)
#	print(res)
#	for i in np.transpose(np.transpose(B)/res):
#		print(i,np.linalg.norm(i))
#	print(np.einsum('ij, ij->i', B,A))
#	coordsubstract(coord1, coord2)
	print (coordsubstract(coord1[0], coord2[0]))
	print (coordsubstract_arr(coord1, coord2))
	sys.exit()
#	distlst = [(23,45,1.45),(2,5,8.5),(9,4,2.2)]
#	distlst.sort(key=lambda tup: tup[2])
#	print distlst
#	array = np.asarray([[1,3,4,2], [6,3,34,5], [6,9,3,1]])
#	print np.ndarray.mean(array, axis=0)
#	print np.ndarray.mean(array, axis=1)
#	
#	satofAx = random.normal(loc=0, scale=.2, size=1000)
#	satofAy = random.normal(loc=0, scale=6.0, size=1000)
#	satofAz = random.normal(loc=0, scale=3.0, size=1000)

#	satofBx = random.normal(loc=0, scale=.2, size=1000)
#	satofBy = random.normal(loc=0, scale=5.0, size=1000)
#	satofBz = random.normal(loc=0, scale=5.0, size=1000)


	satofAx, satofAy, satofAz = random.multivariate_normal([0,0,0], [[0, 2, 5], [2, 4, 0],  [5, 0, 6]], [1000,]).T
	
	satofBx, satofBy, satofBz = random.multivariate_normal([0,0,0], [[1, 5, 1], [5, 4, 1],  [1, 1, 3]], [1000,]).T
	
	for i in range(len(satofBy)):
		satofAx[i] += 7.
			
	sathalos = [satofAx, satofAy, satofAz, satofBx, satofBy, satofBz]
	plotlocal(sathalos, "before2.png")

	for i in range(len(satofBy)):
		satofAx[i] -= 7.
			
	satofAx, satofAy, satofAz = findprince(satofAx, satofAy, satofAz)
	satofBx, satofBy, satofBz = findprince(satofBx, satofBy, satofBz)

	for i in range(len(satofBy)):
		satofAx[i] += 7.
			
	sathalos = [satofAx, satofAy, satofAz, satofBx, satofBy, satofBz]
	plotlocal(sathalos, "after2.png")

#	mean = [0, 0]
#	cov = [[1, 0], [0, 100]]  # diagonal covariance

#	x, y = np.random.multivariate_normal(mean, cov, 5000).T
#	print x.shape, y.shape
#	plt.plot(x, y, 'x')
#	plt.axis('equal')
#	plt.show()

