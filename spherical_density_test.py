from mpl_toolkits.mplot3d import Axes3D
from numpy import sqrt, fabs, loadtxt, arctan2, pi, sin, cos
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter, MaxNLocator, AutoMinorLocator
from constant import *
import os
import sys
from helpers import *
from findprince import *
from histplot import *

def pointinsphere():
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_subplot(111)
	N = 500000
	r = np.random.rand(N)
	theta = np.random.rand(N)*pi
	phi = np.random.rand(N)*pi
	
	x = r * sin(phi) * cos(theta)
	y = r * sin(phi) * sin(theta)
	z = r * cos(phi)
		
	xbins = np.linspace(-1.,1.,200)
	ybins = np.linspace(-1.,1.,200)
	H, xedges, yedges = np.histogram2d(x, z, bins=(ybins, xbins))
	im = ax.imshow(H, extent=[-1.,1.,-1.,1.], interpolation='nearest', aspect = 'equal', vmin=0, vmax=300)
	plt.colorbar(im, format='$%.2f$',cax = fig.add_axes([0.91, 0.2, 0.02, 0.6]))
	plt.savefig('random_sphere.png', transparent=True, dpi=500,bbox_inches='tight')
	plt.close()
	
	var = cos(phi)
#	simplehist(var, 'N(sin phi)', 'NumberOfSatellitePerSinsphi_isotropic_phi.png', 50)
	simplehist(var, 'N(cos phi)', 'NumberOfSatellitePerCosphi_isotropic_phi.png', 50)
#	simplehist(var, 'N(cos theta)', 'NumberOfSatellitePerCostheta_spherical_theta.png', 50)

def pointinsphere2():
	# the generated density is not uniform on the x-z plane because it is a projection onto the plane from 3D distribution
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_subplot(111)
	N = 5000000
	r3d = np.random.rand(N)
#	theta = np.random.rand(N)*pi
	r = np.cbrt(r3d)

	phi = np.random.rand(N)*2e0*pi

	costheta = np.random.rand(N)*2e0-1e0

	x = r * costheta 
	z = r * sqrt(1e0 - costheta ** 2e0) * cos(phi) 
	
	xbins = np.linspace(-1.,1.,200)
	ybins = np.linspace(-1.,1.,200)
	H, xedges, yedges = np.histogram2d(z, x, bins=(ybins, xbins))
	im = ax.imshow(H, extent=[-1.,1.,-1.,1.], interpolation='nearest', aspect = 'equal', vmin=0, vmax=200)
	plt.colorbar(im, format='$%.2f$',cax = fig.add_axes([0.91, 0.2, 0.02, 0.6]))
	plt.savefig('random_sphere.jpg', transparent=True, dpi=500,bbox_inches='tight')
	plt.close()
	
if __name__ == "__main__":
#	pointinsphere()
	pointinsphere2()	
