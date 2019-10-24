import os
import matplotlib.pyplot as plt
from helpers import *
from constant import *
from numpy import loadtxt, sqrt
import sys
'''
plotting the positions and mass of the small satellites near the control (artificially Local Group binaries imitating the real binaries)
'''

def plotlocalprince(pairno, mA, mB, hostApos, hostBpos, sathalos, satofAm, satofBm, princeaxesA, princeaxesB, plotname):
	mid = middlecoordmpc(hostApos, hostBpos)
	dsep = coorddistance(hostApos, hostBpos)
	plt.clf()
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_subplot(221, projection='3d')

	[satofAx, satofAy, satofAz, satofBx, satofBy, satofBz] = sathalos
	
	for i in range(len(satofAx)):
		mass = satofAm[i]
		ax.scatter(satofAx[i], satofAy[i], satofAz[i], zdir='z', color = 'black', s = mass/plotmassscale)

	ax.set_xlim3d(mid[0]-dsep, mid[0]+dsep)
	ax.set_ylim3d(mid[1]-dsep, mid[1]+dsep)	
	ax.set_zlim3d(mid[2]-dsep, mid[2]+dsep)
	
	for i in range(len(satofBx)):
		mass = satofBm[i]
		ax.scatter(satofBx[i], satofBy[i], satofBz[i], zdir='z', color = 'black', s = mass/plotmassscale)

	shiftedaxesA = [[],[],[]]
	shiftedaxesB = [[],[],[]]
	for i in range(3):
		for j in range(3):
			shiftedaxesA[i].append(princeaxesA[i][j]/5 + hostApos[j]) #5 for scaling length of the eigenvector
#			print princeaxesA[i][j], hostApos[i]
			shiftedaxesB[i].append(princeaxesB[i][j]/5 + hostBpos[j])

	for i in range(3):
		ax.plot([hostApos[0], shiftedaxesA[i][0]], [hostApos[1], shiftedaxesA[i][1]], [hostApos[2], shiftedaxesA[i][2]], color = 'purple')
		ax.plot([hostBpos[0], shiftedaxesB[i][0]], [hostBpos[1], shiftedaxesB[i][1]], [hostBpos[2], shiftedaxesB[i][2]], color = 'green')
		
	#take care of periodic boundary condition for binary
	hostApos, hostBpos, marker = coordwrap(hostApos, hostBpos) #if marker[i] = 1, hostApos 63.5 -> -.5; marker[i] = 2 hostBpos 63.5 -> -.5
#	ax.scatter(hostApos[0], hostApos[1], hostApos[2], alpha = 0.3, color = 'red', s = mA/plotmassscale/10) 	#halo size proportional to mass
#	ax.scatter(hostBpos[0], hostBpos[1], hostBpos[2], alpha = 0.3, color = 'blue', s = mB/plotmassscale/10)
#	
	ax2 = fig.add_subplot(222)

	ax2.set_xlim(mid[0]-dsep, mid[0] +dsep)
	ax2.set_ylim(mid[1]-dsep, mid[1] +dsep)	
	
	for i in range(len(satofAx)):
		mass = satofAm[i]
		ax2.scatter(satofAx[i], satofAy[i], color = 'black', s = mass/plotmassscale)

	for i in range(len(satofBx)):
		mass = satofBm[i]
		ax2.scatter(satofBx[i], satofBy[i], color = 'black', s = mass/plotmassscale)
#	ax2.scatter(hostApos[0], hostApos[1], color = 'red', s = mA/plotmassscale/10) 	#halo size proportional to mass
#	ax2.scatter(hostBpos[0], hostBpos[1], color = 'blue', s = mB/plotmassscale/10)

	colorlist = ['purple', 'green', 'blue']

	for i in range(3):
		ax2.plot([hostApos[0], shiftedaxesA[i][0]], [hostApos[1], shiftedaxesA[i][1]], color = colorlist[i])
		ax2.plot([hostBpos[0], shiftedaxesB[i][0]], [hostBpos[1], shiftedaxesB[i][1]], color = colorlist[i])
		
	ax3 = fig.add_subplot(223)
	ax3.set_xlim(mid[1]-dsep, mid[1] +dsep)	
	ax3.set_ylim(mid[2]-dsep, mid[2] +dsep)
	
	for i in range(len(satofAx)):
		mass = satofAm[i]
		ax3.scatter(satofAy[i], satofAz[i], color = 'black', s = mass/plotmassscale)

	for i in range(len(satofBx)):
		mass = satofBm[i]
		ax3.scatter(satofBy[i], satofBz[i], color = 'black', s = mass/plotmassscale)
#	ax3.scatter(hostApos[1], hostApos[2], color = 'red', s = mA/plotmassscale/10) 	#halo size proportional to mass
#	ax3.scatter(hostBpos[1], hostBpos[2], color = 'blue', s = mB/plotmassscale/10)
	for i in range(3):
		ax3.plot([hostApos[1], shiftedaxesA[i][1]], [hostApos[2], shiftedaxesA[i][2]], color = colorlist[i])
		ax3.plot([hostBpos[1], shiftedaxesB[i][1]], [hostBpos[2], shiftedaxesB[i][2]], color = colorlist[i])	
	ax4 = fig.add_subplot(224)

	ax4.set_xlim(mid[0]-dsep, mid[0] +dsep)
	ax4.set_ylim(mid[2]-dsep, mid[2] +dsep)
	
	for i in range(len(satofAx)):
		mass = satofAm[i]
		ax4.scatter(satofAx[i], satofAz[i], color = 'black', s = mass/plotmassscale)

	for i in range(len(satofBx)):
		mass = satofBm[i]
		ax4.scatter(satofBx[i], satofBz[i], color = 'black', s = mass/plotmassscale)
#	ax4.scatter(hostApos[0], hostApos[2], color = 'red', s = mA/plotmassscale/10) 	#halo size proportional to mass
#	ax4.scatter(hostBpos[0], hostBpos[2], color = 'blue', s = mB/plotmassscale/10)		
	for i in range(3):
		ax4.plot([hostApos[0], shiftedaxesA[i][0]], [hostApos[2], shiftedaxesA[i][2]], color = colorlist[i])
		ax4.plot([hostBpos[0], shiftedaxesB[i][0]], [hostBpos[2], shiftedaxesB[i][2]], color = colorlist[i])	
	fig.savefig(plotname, dpi=500)
	plt.close()
		
