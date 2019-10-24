from scipy.stats import chisquare, ks_2samp
import os
from numpy import loadtxt, arctan2, sqrt, cos, fabs, sqrt
import matplotlib.pyplot as plt
import numpy as np
import sys
from density2dplot import *
from constant import *
from histplot import *

def Nplot_relDistSat(r, rf, rf2):
	fig = plt.figure(figsize=(6,6))
	ax = fig.add_subplot(111)

	ax.set_xlabel("$d_{\mathrm{sat}}/d_{\mathrm{sep}}$",fontsize=14)
	ax.set_ylabel("$N_{\mathrm{sat}}$",fontsize=14)
	n_bins = 50
	plt.hist(r, n_bins, histtype = 'step', cumulative=True, label = 'signal')
	plt.hist(rf, n_bins, histtype = 'step', cumulative=True, label = 'control w/ $R_{[\mathrm{ol}]}$ = 1.5$d_{\mathrm{sep}}$')	
	plt.hist(rf2, n_bins, histtype = 'step', cumulative=True, label = 'control w/ $R_{[\mathrm{ol}]}$ = 1$d_{\mathrm{sep}}$')	
	
	legend = ax.legend(loc="upper center", shadow=True, fontsize = 8)
	fig.savefig("Ndist_relDsat_signvsOL.pdf", dpi=500)
	plt.close()

def twoDpositionplot(costheta, dist, hostmass, dsep, d):
	print(costheta[:10], dist[:10], hostmass[:10])
	fig = plt.figure(figsize=(6,6))
	ax = fig.add_subplot(111)
	
	B = np.where(hostmass < 1e0)
	A = np.where(hostmass > 1e0)
	Acostheta = costheta[A[0]]
	Bcostheta = costheta[B[0]]
	Adist = dist[A[0]] * dsep[A[0]]
	Bdist = dist[B[0]] * dsep[B[0]]
	
	theta_A = np.random.rand(len(A[0]))*2e0*pi

	theta_B = np.random.rand(len(B[0]))*2e0*pi
	
	x_satA = Adist * Acostheta
	y_satA = Adist * sqrt(1e0-Acostheta**2e0) * cos(theta_A)

	x_satB = dsep[B[0]] - Bdist * Bcostheta
	y_satB = Bdist * sqrt(1e0-Bcostheta**2e0) * cos(theta_B)
	
	x_sat = np.concatenate((x_satA, x_satB))
	y_sat = np.concatenate((y_satA, y_satB))
	print(x_sat[:10])
	print(y_sat[:10])	
#	xmax = 4e0
#	xmin = -4e0
#	ymax = 4e0
#	ymin = -4e0
#					
	xmax = 2e0
	xmin = -1e0
	ymax = 1e0
	ymin = -1e0
	nbins = 40
	vmin0 = 1.6
	vmax0 = 4
		
	deltax = (xmax-xmin)/(2e0 * nbins)
	deltay = (ymax-ymin)/(2e0 * nbins)
	
	x0 = np.arange(xmin, xmax + deltax, deltax)
	y0 = np.arange(ymax, ymin - deltay, -deltay)
		
	x = np.arange(xmin + deltax/2e0, xmax, deltax) # plot grids center
	y = np.arange(ymax - deltay/2e0, ymin, -deltay)
	X, Y = np.meshgrid(x,y)
	
	density_mtrA = []
	for i in range(nbins*2): #x direction grid
		density_mtrA.append([])
		for k in range(nbins*2): #y direction grid
			index = np.where( np.logical_and(np.logical_and(x_satA > x0[k], x_satA < x0[k + 1]), np.logical_and(y_satA < y0[i], y_satA > y0[i + 1])))
#			print(x0[k], x0[k+1],index[0])
#			print(y0[k], y0[k+1],index[0])
			density_mtrA[i].append(np.log10(len(index[0])))

	im = ax.imshow(np.array(density_mtrA), extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Reds, aspect = 'equal', interpolation='nearest', vmin=vmin0, vmax=vmax0, alpha = 0.7)
	ax.contour(X, Y, np.array(density_mtrA), colors='red', alpha = 0.3,levels =np.arange(.6,4.,.4))
	
	density_mtrB = []
	for i in range(nbins*2): #x direction grid
		density_mtrB.append([])
		for k in range(nbins*2): #y direction grid
			index = np.where( np.logical_and(np.logical_and(x_satB > x0[k], x_satB < x0[k + 1]), np.logical_and(y_satB < y0[i], y_satB > y0[i + 1])))
			density_mtrB[i].append(np.log10(len(index[0])))

	im = ax.imshow(np.array(density_mtrB), extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Blues, aspect = 'equal', interpolation='nearest', vmin=vmin0, vmax=vmax0, alpha = 0.7)
	ax.contour(X, Y, np.array(density_mtrB), colors='blue', alpha = 0.3,levels =np.arange(.6,4.,.4))

#	density_mtr = []
#	for i in range(nbins*2): #x direction grid
#		density_mtr.append([])
#		for k in range(nbins*2): #y direction grid
#			index = np.where( np.logical_and(np.logical_and(x_sat > x0[k], x_sat < x0[k + 1]), np.logical_and(y_sat < y0[i], y_sat > y0[i + 1])))
#			density_mtr[i].append(np.log10(len(index[0])))

##		im = ax.imshow(np.array(density_mtr), extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Blues, aspect = 'equal', interpolation='nearest', vmin=vmin0, vmax=vmax0)
#	im = ax.imshow(np.array(density_mtr), extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Blues, aspect = 'equal', interpolation='nearest')		
#	ax.contour(X, Y, np.array(density_mtr), colors='black')

	plt.colorbar(im, format='$%.2f$', ax = ax)
	ax.grid(True, linestyle= '-', which= 'major',color= '0.75')
	ax.grid(True, which='both')
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	
#	t = np.arange(vmin0,vmax0,.4)
#	fig.colorbar(im, format='$%.2f$',cax = fig.add_axes([0.91, 0.2, 0.02, 0.6]))
	fig.savefig('ABdensity_control_2color_zoomout_dsep_%.1f_to_%.1f.pdf'% (d, d+0.1), dpi=500)			
	plt.close()

def plotNvsmassratio(HostMassRatio):
	nbins = 40
	
	MRmax = 1.0
	Massratiobins = np.arange(0, MRmax + MRmax/nbins, MRmax/nbins)

	for i in range(len(HostMassRatio)):
		if HostMassRatio[i] >1:
			HostMassRatio[i] = 1e0/HostMassRatio[i]
					
	histx, bin_edgesx = np.histogram(HostMassRatio, bins=Massratiobins, normed = True)
	
	MRlocs = np.arange(0 + .5 * MRmax/nbins, MRmax + .5*MRmax/nbins, MRmax/nbins)
			
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')	
	fig = plt.figure()
	fig.set_size_inches(8, 5, forward=True)
	ax = fig.add_subplot(111)
			
	ax.plot(MRlocs, histx, color='blue', linewidth = 1.5)
	ax.set_xlim(0, MRmax)	
	ax.set_xlabel('$M_1/M_2$', fontsize=20)
	ax.set_ylabel('$N(>M_1/M_2)$',fontsize=20)
	plt.grid(True)
	fig.savefig('N_vs_Massratio_fake.jpg', dpi=500, bbox_inches='tight')

def plotcontrol_atlowmassratio(cos_fake, Hostmassratio_fake):
	nbins = 40

	for i in range(len(Hostmassratio_fake)):
		if Hostmassratio_fake[i] > 1:
			Hostmassratio_fake[i] = 1e0/Hostmassratio_fake[i]

	smallHostMassRatio = np.where(Hostmassratio_fake < 0.1)
	cos = cos_fake[smallHostMassRatio[0]]
	cosbins = np.arange(-1.0, 1.0 + (1.0-(-1.0))/nbins, (1.0-(-1.0))/nbins)
	coslocs = np.arange(-1.0 + .5 * (1.0-(-1.0))/nbins, 1.0 + .5*(1.0-(-1.0))/nbins, (1.0-(-1.0))/nbins)
		
	histx, bin_edgesx = np.histogram(cos, bins=cosbins, normed = True)
	print(histx)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')	
	fig = plt.figure()
	fig.set_size_inches(8, 5, forward=True)
	ax = fig.add_subplot(111)
	
	indpos = np.where(np.array(cos) > 0.8)
	indneg = np.where(np.array(cos) < -0.8)	
	
	Qplus90 = len(indpos[0])
	Qminus10 = len(indneg[0])

	meanAB, sigmaAB = markovN(len(cos))
	meanabs, sigmaabs = markovNabs(len(cos))
	
	ax.text(-.5, 0.66, "$q=%.3f, %.0f \sigma$" % (Qplus90/Qminus10, abs((Qplus90/Qminus10 - meanAB)/sigmaAB)), color='blue', fontsize=18, fontweight='bold')			
	ax.plot(coslocs, histx, color='blue', linewidth = 1.5)
	ax.set_xlim(-1, 1)	
	ax.set_ylim(0, 0.6)	
	ax.set_xlabel('$\\cos\\theta$', fontsize=20)
	ax.set_ylabel('$p(\\cos\\theta)$ for $M2/M1 < 0.1$',fontsize=20)
	plt.grid(True)
	fig.savefig('smallMassratio_control.jpg', dpi=500, bbox_inches='tight')

def plotcontrol_athighmassratio(cos_fake, Hostmassratio_fake):
	nbins = 40

	for i in range(len(Hostmassratio_fake)):
		if Hostmassratio_fake[i] > 1:
			Hostmassratio_fake[i] = 1e0/Hostmassratio_fake[i]

	smallHostMassRatio = np.where(np.logical_and(Hostmassratio_fake > 0.3, Hostmassratio_fake < 0.4))
	cos = cos_fake[smallHostMassRatio[0]]
	cosbins = np.arange(-1.0, 1.0 + (1.0-(-1.0))/nbins, (1.0-(-1.0))/nbins)
	coslocs = np.arange(-1.0 + .5 * (1.0-(-1.0))/nbins, 1.0 + .5*(1.0-(-1.0))/nbins, (1.0-(-1.0))/nbins)
		
	histx, bin_edgesx = np.histogram(cos, bins=cosbins, normed = True)
	print(histx)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')	
	fig = plt.figure()
	fig.set_size_inches(8, 5, forward=True)
	ax = fig.add_subplot(111)

	indpos = np.where(np.array(cos) > 0.8)
	indneg = np.where(np.array(cos) < -0.8)	
	
	Qplus90 = len(indpos[0])
	Qminus10 = len(indneg[0])

	meanAB, sigmaAB = markovN(len(cos))
	meanabs, sigmaabs = markovNabs(len(cos))
	
	ax.text(-.5, 0.66, "$q=%.3f, %.0f \sigma$" % (Qplus90/Qminus10, abs((Qplus90/Qminus10 - meanAB)/sigmaAB)), color='blue', fontsize=18, fontweight='bold')	
				
	ax.plot(coslocs, histx, color='blue', linewidth = 1.5)
	ax.set_xlim(-1, 1)	
	ax.set_ylim(0, 0.6)	
	ax.set_xlabel('$\\cos\\theta$', fontsize=20)
	ax.set_ylabel('$p(\\cos\\theta)$ for $M2/M1 < 0.1$',fontsize=20)
	plt.grid(True)
	fig.savefig('largeMassratio_control.jpg', dpi=500, bbox_inches='tight')

def chisqrtest(coslist, coslist_fake):
	for binnum in [40,80,100,200, 500, 1000]:
		histx1, bin_edgesx1 = np.histogram(coslist, bins=binnum, range = (-1,1))
		histx2, bin_edgesx2 = np.histogram(coslist_fake, bins=binnum, range = (-1,1))
		histx1 = histx1/sum(histx1)
		histx2 = histx2/sum(histx2)
		
		K1 = np.sqrt(sum(histx2)/sum(histx1))
		K2 = np.sqrt(sum(histx1)/sum(histx2))
#		print(K1, K2)
		chisqr = 0
		for i in range(len(histx1)):
			chisqr += (K1 * histx1[i] - K2 * histx2[i])**2e0/(histx1[i] + histx2[i])
			
#		print(histx1, histx2)
		print(binnum, chisqr)
#	print(ks_2samp(coslist, coslist_fake))
#	print(coslist[:100])
#	print(coslist_fake[:100])
	
if __name__ == "__main__":
	#loading signal
	x1,y1 = loadtxt(hrcosine_distn_file, skiprows=1, usecols= (0, 10), unpack=True)
#	simplehist(x1, "cos", "cos2.png", 40)
	#loading overlapping artificial signal
	x2,y2 = loadtxt(hrcosine_distn_file_fake, skiprows=1, usecols= (0,7), unpack=True)
	
	mass = loadtxt(hrcosine_distn_file, skiprows=1, usecols= (6,), unpack=True)
	mass2 = loadtxt(hrcosine_distn_file_fake, skiprows=1, usecols= (8,), unpack=True)
	Hostmassratio = loadtxt(hrcosine_distn_file_fake, skiprows=1, usecols= (1,), unpack=True)
	
	control_dsep = loadtxt(hrcosine_distn_file_fake, skiprows=1, usecols= (5,), unpack=True)
	
	dsep_grid = np.arange(0.3,1.5,0.1)
	
	toosmallhalo1 = (mass /particlemass < 20)
	toosmallhalo2 = (mass2 /particlemass < 20)	
	idofbigenoughhalos1 = np.where(toosmallhalo1 == 0)
	idofbigenoughhalos2 = np.where(toosmallhalo2 == 0)
	
	x1 = x1[idofbigenoughhalos1[0]]
	y1 = y1[idofbigenoughhalos1[0]]
	x2 = x2[idofbigenoughhalos2[0]]
	y2 = y2[idofbigenoughhalos2[0]]
	control_dsep = control_dsep[idofbigenoughhalos2[0]]
	Hostmassratio_fake = Hostmassratio[idofbigenoughhalos2[0]]
	
#	plotNvsmassratio(Hostmassratio_fake)
	
#	plotcontrol_atlowmassratio(x2, Hostmassratio_fake)
#	plotcontrol_athighmassratio(x2, Hostmassratio_fake)
#	sys.exit()

##	id_smallHMR = np.where( Hostmassratio > 0.8 or Hostmassratio < 1.2)
##	id_unityHMR = np.where( np.abs(Hostmassratio-0.2) <1.)
#	id_largeHMR = np.where( Hostmassratio > 1.2 )

#	id_small_dsep = np.where( control_dsep < 0.4)
#	id_large_dsep = np.where( control_dsep > 0.6)
	
#	print(len(x2), len(x2[idofbigenoughhalos2]))
#	sys.exit()

#	chisqrtest(x1, x2)
#	sys.exit()
	
	plotname = "rtheta_densityplot_signal_3panels_40bins.pdf"
	densityplot_2panels_signal(x1, y1, x2, y2, '$\cos \\theta$', '$d_{\mathrm{sat}}/d_{\mathrm{sep}}$', os.path.join(os.getcwd(), "paper", plotname ), [-1e0,1e0], [0,.5e0], [0.4, .7], 2.5, 40)
	print(len(x1), len(x2))
	sys.exit()
	
	for i, d in enumerate(dsep_grid[:-1]):
		id0 = np.where( np.logical_and(control_dsep > d, control_dsep < dsep_grid[i+1]))
		plotname = "rtheta_densityplot_signal_3panels_40bins_dsep_%.1f_to_%.1f_samplesize%.2f.pdf" % (d, d+0.1, len(id0[0])/len(x1))
		densityplot_3panels_signal(x1[id0[0]], y1[id0[0]], x2[id0[0]], y2[id0[0]], '$\cos \\theta$', '$d_{\mathrm{sat}}/d_{\mathrm{sep}}$', os.path.join(os.getcwd(), "paper", plotname ), [-1e0,1e0], [0,.5e0], [0.4, .7], 2.3, 40)
#		twoDpositionplot(x2[id0[0]], y2[id0[0]], Hostmassratio[id0[0]], control_dsep[id0[0]], d)
#	densityratioplot_2panels(x1, y1, x2, y2, x3, y3, '$\cos \\theta$', '$d_{\mathrm{sat}}/d_{\mathrm{sep}}$', os.path.join(os.getcwd(), "paper", "rtheta_densityplot_ratioplot_40bins.pdf"), [-1e0,1e0], 40)	
