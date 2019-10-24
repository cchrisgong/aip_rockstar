'''Roche potential fitted to satellite probability density of local group pair'''
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from numpy import sqrt, fabs, loadtxt, arctan2, pi, log, amin
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter, MaxNLocator, AutoMinorLocator
from constant import *
import os
import sys
from helpers import *
from findprince import *
from histplot import *
from matplotlib import rcParams
import matplotlib as mpl

def lawofcos(coordS, coordA, coordB):
	b_sq = coorddistance_sqd(coordS, coordA)
	c_sq = coorddistance_sqd(coordB, coordA)
	a_sq = coorddistance_sqd(coordS, coordB) #B needs to be the host further to S (satellite) in this function
	if b_sq == 0 or c_sq == 0:
		return "small satellite distance to one of the hosts is rounded to zero"
	else:
		return (b_sq + c_sq - a_sq)/(2e0* sqrt(b_sq)*sqrt(c_sq))
		
def rochePot(q, x, y, z):
	#x y z in ratio of true x y z of satelllite with dsep
	r1 = sqrt(( x - 1e0 ) ** 2e0 + y ** 2e0 + z ** 2e0)
	r2 = sqrt(( x ) ** 2e0 + y ** 2e0 + z ** 2e0)
	qfac = q / (1e0 + q)
	return 2e0 * q / ((1e0 + q) * r1) + 2e0 / ((1e0 + q) * r2) + (x - qfac) ** 2e0 + y ** 2e0
	
def drawRochePot():
	xplotrange = np.linspace(-.5e0, 1.5e0, 300)
	yplotrange = np.linspace(-1e0, 1e0, 300)	
	z = 0
	fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8)) = plt.subplots(2, 4, sharex='col', sharey='row')
	fig.set_size_inches(11, 5, forward = True)		
	ax_list = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]	

	for i in range(len(mlist)):
		x, y = np.meshgrid(xplotrange, yplotrange)
		cs = ax_list[i].contour(x, y, rochePot(mlist[i], x, y, z), levels= np.logspace(0e0,1e0,50))
		ax_list[i].set_xlabel("$M_{2}/M_{1}$ = " + str(mlist[i]))
		ax_list[i].xaxis.set_major_locator(MaxNLocator(3))	
		
	plt.colorbar(cs, format='$%.2f$',cax = fig.add_axes([0.91, 0.15, 0.02, 0.7]))
	plt.savefig("roche2.png", dpi=500, bbox_inches='tight')

def cos_satellite_angle_rotated(mA, mB, hostApos, hostBpos, dsep, RvirA, RvirB, SMALL, j):
	################## prepare satellite list ###############
	[SMALLm, SMALLx, SMALLy, SMALLz] = SMALL
	
	small = []
	for i in range(len(SMALLm)):
#		if SMALLm[i] < min(mA, mB): #make sure the hosts are not counted as its own satellite
		small.append([SMALLm[i], SMALLx[i], SMALLy[i], SMALLz[i]]) 

	#using half distance
	rA, rB = dsep/2e0, dsep/2e0

	satofA, satofB = [], []
	satofAm, satofBm = [], []

	for i in range(len(small)):
		pos = np.array(small[i][1:])
		Sm = np.array(small[i][0])

		satDistA = coorddistance(pos, hostApos)
		satDistB = coorddistance(pos, hostBpos)

		if satDistA < rA: #sat is closer to host A, belongs to host A
			coorddiff = coordsubstract(pos, hostApos)
			satofA.append(coorddiff)
#			satofA.append(pos)
			satofAm.append(Sm)
		
		if satDistB < rB: #sat is closer to host B, belongs to host B
			coorddiff = coordsubstract(pos, hostBpos)
			satofB.append(coorddiff)
#			satofB.append(pos)		
			satofBm.append(Sm)
	
	if coordcompare(hostBpos[0], hostApos[0]): #B is to the right of A in a wrap box in x direction
		vAB = coordsubstract(hostBpos,hostApos)
		tag = 1

	if coordcompare(hostBpos[0], hostApos[0]) == 0:
		vAB = coordsubstract(hostApos, hostBpos)
		tag = 0

	vAB = coordsubstract(hostApos, hostBpos)
	vAB_norm = np.linalg.norm(vAB)
#	print(vAB)
	vAB = vAB/vAB_norm
	xaxis1 = [1e0, 0e0, 0e0]
	costheta1 = np.dot(vAB, xaxis1)
	sintheta1 = np.linalg.norm(np.cross(vAB, xaxis1))
	theta1 = arctan2(sintheta1,costheta1)

	xaxis2 = [-1e0, 0e0, 0e0]
	costheta2 = np.dot(vAB, xaxis2)
	sintheta2 = np.linalg.norm(np.cross(vAB, xaxis2))
	theta2 = arctan2(sintheta2,costheta2)

	if theta1 > pi/2e0:
		xaxis = xaxis2
		theta = theta2
		costheta = costheta2
		sintheta = sintheta2
		
	if theta2 > pi/2e0:
		xaxis = xaxis1
		theta = theta1
		costheta = costheta1
		sintheta = sintheta1

	if tag == 1:
		[v1, v2, v3] = np.cross(xaxis,vAB)
		Anewpos = np.array([0, 0, 0])
		Bnewpos = np.array([dsep, 0, 0])
		
	if tag == 0:
		[v1, v2, v3] = np.cross(xaxis, vAB)
		Anewpos = np.array([dsep, 0, 0])
		Bnewpos = np.array([0, 0, 0])
		
	vmat = np.asarray([[0e0, -v3, v2], [v3, 0e0, -v1], [-v2, v1, 0e0]])
	vmatXvmat = np.asarray([[v1*v1, v1*v2, v1*v3], [v1*v2, v2*v2, v2*v3], [v1*v3, v2*v3, v3*v3]])
	RotMatrix = costheta * np.identity(3) + sintheta * vmat + vmatXvmat * (1e0 - costheta)
	
	for i in range(len(satofA)):
		satofA[i] = np.dot(satofA[i], RotMatrix)
	for i in range(len(satofB)):
		satofB[i] = np.dot(satofB[i], RotMatrix)
	#########################################################
	
	cosA = []
	cosB = []
	
	for i in range(len(satofA)):
		if tag == 0:
			satofA[i][0] += dsep
		pos2 = np.array(satofA[i])
#		if fabs(pos2[2]) > 0. * dsep and fabs(pos2[2]) < 0.05 * dsep:
		satDistA = coorddistance(pos2, Anewpos)
		Sm = satofAm[i]
		cosangle = lawofcos(pos2, Anewpos, Bnewpos)
		if isinstance(cosangle, str) == 0:
			cosA.append([cosangle, Sm, satDistA, pos2])
	
	for i in range(len(satofB)):
		if tag == 1:
			satofB[i][0] += dsep
		pos2 = np.array(satofB[i])
#		if fabs(pos2[2]) > 0. * dsep and fabs(pos2[2]) < 0.05 * dsep:
		satDistB = coorddistance(pos2, Bnewpos)
		Sm = satofBm[i]
		cosangle = lawofcos(pos2, Bnewpos, Anewpos)
		if isinstance(cosangle, str) == 0:
			cosB.append([cosangle, Sm, satDistB, pos2])
	
	return cosA, cosB, satofA, satofB, satofAm, satofBm, RotMatrix, Anewpos, Bnewpos, tag

def cosdistn(LGcathaloID, LGm, LGx, LGy, LGz, LGRvir, LGb, LGc):
	#writing to file the small satellite and their host properties
	
	f = open(hrcosine_distn_file_roche,"w")
##	 cos of angle from sat to Host		host/partner_massratio		eccentricity c/a(for binary) 	eccentricity b/a(for binary)	eccentricity c/b (for binary)	mass of satellite/mass of host		mass of satellite/mass of partner 		mass of satellite 		distance of binary 		distance from satellite to host		distance from satellite to host/distance of binary
	f.write("hostID" +  "\t" + "cos_angle" +  "\t" + "host/partner_massratio" + "\t" + "c/a" + "\t" + "b/a" + "\t" + "c/b" + "\t" + "mS/mH" + "\t" + "mS/mP" + "\t" + "mS" + "\t" + "dH"  + "\t" + "dS" + "\t" + "dS/dH" + "\n")
	neglist = loadtxt("negfile.dat")
	for j in range(len(LGm)):
		if j % 2 == 0:
			pairno = int(j/2e0) #unpack the data stored in "low_res_trueLG_pop.dat" in pairs (pairs were written sequentially)
			
			if pairno not in neglist:
				small_halo_path = os.path.join(os.getcwd(), folder_smallhalo, str(pairno) + "_small_halos_feats.dat")
				SMALLm, SMALLx, SMALLy, SMALLz = loadtxt(small_halo_path, skiprows = 1, usecols = (masscol, xcol, ycol, zcol), unpack=True)#get small halo positions
				
				SMALL = [SMALLm, SMALLx, SMALLy, SMALLz]
				
#				mA = bigger partner
#				mB = smaller partner

				mA = max([LGm[j], LGm[j+1]])
				mB = min([LGm[j], LGm[j+1]])
				indA = j+[LGm[j], LGm[j+1]].index(mA)
				indB = j+[LGm[j], LGm[j+1]].index(mB)
				bA, bB = LGb[indA], LGb[indB] #host eccentricity; b is b/a
				cA, cB = LGc[indA], LGc[indB] #host eccentricity; c is c/a
				hostApos = np.array([LGx[indA], LGy[indA], LGz[indA]])
				hostBpos = np.array([LGx[indB], LGy[indB], LGz[indB]])
				AcatID = LGcathaloID[indA]
				BcatID = LGcathaloID[indB]
				dsep = coorddistance(hostApos, hostBpos)
				RvirA, RvirB = LGRvir[indA]/1000e0, LGRvir[indB]/1000e0 #host radius

				print ("pair number", pairno)
				cosA, cosB, satofA, satofB, satofAm, satofBm, princeaxesAB, hostApos, hostBpos, ABtag = cos_satellite_angle_rotated(mA, mB, hostApos, hostBpos, dsep, RvirA, RvirB, SMALL, j) #ABtag: 1, B is on the right
				
#				if j > 60 and j < 100:
#					vAB1 = coordsubstract(hostBpos, hostApos)
#					sathalos = [satofA, satofB]
#					plotname = "HRpair" + str(pairno) + "_wsmallhalo_ABconnected_unshifted.png"
##					plotname = "HRpair" + str(pairno) + "_wsmallhalo_ABconnected_before_unshifted.png"
#					plotlocalprince(pairno, mA, mB, sathalos, satofAm, satofBm, princeaxesAB, vAB1, hostApos, hostBpos, plotname, dsep, 1e0)
#				if j > 100:
#					sys.exit()
				
#	#			note: cosA[i] = [cosangle, mass of satellite, sat dist to host, position of satellite]
				for i in range(len(cosA)):
					f.write(str(AcatID) + "\t" + str(cosA[i][0]) + "\t" + str(mA/mB) + "\t" + str(cA) + "\t" + str(bA) + "\t" + str(cA/bA) + "\t" + str(cosA[i][1]/mA) + "\t" + str(cosA[i][1]/mB) + "\t" + str(cosA[i][1]) + "\t" + str(cosA[i][3][0])+ "\t" + str(cosA[i][3][1])+ "\t" + str(cosA[i][3][2]) + "\t" + str(dsep) + "\t" + str(cosA[i][2]) + "\t" + str(cosA[i][2]/dsep) +"\t" + str(ABtag) + "\n")

				for i in range(len(cosB)):
					f.write(str(BcatID) + "\t" + str(cosB[i][0]) + "\t" + str(mB/mA) + "\t" + str(cB) + "\t" + str(bB) + "\t" + str(cB/bB) + "\t" + str(cosB[i][1]/mB) + "\t" + str(cosB[i][1]/mA) + "\t" + str(cosB[i][1]) + "\t" + str(cosB[i][3][0])+ "\t" + str(cosB[i][3][1])+ "\t" + str(cosB[i][3][2]) + "\t" + str(dsep) + "\t" + str(cosB[i][2]) + "\t" + str(cosB[i][2]/dsep) + "\t" + str(ABtag) + "\n")
	f.close()

def drawdata_onepanel():
	hostcatID, coslist, HostMassRatio, covera, bovera, coverb, sToHmass, sToPmass, sMass, sPosx, sPosy, sPosz, dsep, distsat, relDistSat, tag = loadtxt(hrcosine_distn_file_roche, skiprows=1, unpack=True)
	
	toosmallhalo = (sMass /particlemass < 20)		
	idofbigenoughhalos = np.where(toosmallhalo == 0)
	
	coslist = coslist[idofbigenoughhalos[0]]
	relDistSat = relDistSat[idofbigenoughhalos[0]]
	dsep = dsep[idofbigenoughhalos[0]]
	HostMassRatio = HostMassRatio[idofbigenoughhalos[0]]
	sPosx = sPosx[idofbigenoughhalos[0]]
	sPosy = sPosy[idofbigenoughhalos[0]]
	sPosz = sPosz[idofbigenoughhalos[0]]
	sMass = sMass[idofbigenoughhalos[0]]
	
	xlistplot = []
	ylistplot = []
	coslistplot =[]
	tol = .05e0
	
	for j in range(len(mlist)):
		xlistplot.append([])
		ylistplot.append([])
		coslistplot.append([])
		for i in range(len(HostMassRatio)):
			if fabs(HostMassRatio[i] - mlist[j]) < tol: # HostMassRatio = host/partner
				if tag[i] == 0: #A
					xlistplot[j].append((-sPosx[i] + dsep[i])/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])
				else: #B
					xlistplot[j].append(sPosx[i]/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])					
				coslistplot[j].append(coslist[i]) #only minor host distn is plotted in histogram
		
			if fabs(1e0/HostMassRatio[i] - mlist[j]) < tol:
				if tag[i] == 0:
					xlistplot[j].append((-sPosx[i] + dsep[i])/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])
				else:
					xlistplot[j].append(sPosx[i]/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])
#		print(len(coslistplot[j]))		

	binno = 40
	plotpath = os.path.join(os.getcwd(),"roche_1panel.png")	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')	
	fig = plt.figure()
	
	fig.set_size_inches(5, 4, forward=True)
	ax = fig.add_subplot(111)
		
	xmin, xmax, ymin,ymax = -.5e0, 1.5e0, -.6e0, .6e0
	xlims, ylims = [xmin,xmax], [ymin,ymax]
	nxbins, nybins, nbins = binno, binno, binno
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	
	cosbins = np.arange(-1.0, 1.0 + (1.0-(-1.0))/nbins, (1.0-(-1.0))/nbins)
	coslocs = np.arange(-1.0 + .5 * (1.0-(-1.0))/nbins, 1.0 + .5*(1.0-(-1.0))/nbins, (1.0-(-1.0))/nbins)
	for i in range(len(mlist)):
		ax.set_xlim(-1, 1)
		ax.set_ylim(0.38, 0.86)
		ticklabels = ax.get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(12)
			
		#Plot the histograms
		histx, bin_edgesx = np.histogram(coslistplot[i], bins=cosbins, normed = True)
		Qminus10 = sum(histx[:int(binno*0.1)])
		Qplus90 = sum(histx[int(binno*0.9):])
		
		meanAB, sigmaAB = markovN(len(coslistplot[i]))
		q = Qplus90/Qminus10
		
		ax.plot(coslocs, histx, linewidth = 1.5, label = " $%.1f $ < M_{2}/M_{1} <$ %.1f, q = %.3f, %.0f \sigma$ " % (mlist[i]-0.05, mlist[i]+0.05, q, (q - meanAB) / sigmaAB))
		
		
		#Make the tickmarks pretty

		ticklabels = ax.get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(15)
				
		# Remove the inner axes numbers of the histograms
		
	ax.set_xlabel("$\cos \\theta$",fontsize=15)
	ax.set_ylabel("$p(\cos \\theta)$",fontsize=15)

	plt.savefig(plotpath, transparent=True, dpi=500,bbox_inches='tight')
	plt.close()
	
def drawdata_onerow():
	hostcatID, coslist, HostMassRatio, covera, bovera, coverb, sToHmass, sToPmass, sMass, sPosx, sPosy, sPosz, dsep, distsat, relDistSat, tag = loadtxt(hrcosine_distn_file_roche, skiprows=1, unpack=True)
	
	toosmallhalo = (sMass /particlemass < 20)		
	idofbigenoughhalos = np.where(toosmallhalo == 0)
	
	coslist = coslist[idofbigenoughhalos[0]]
	relDistSat = relDistSat[idofbigenoughhalos[0]]
	dsep = dsep[idofbigenoughhalos[0]]
	HostMassRatio = HostMassRatio[idofbigenoughhalos[0]]
	sPosx = sPosx[idofbigenoughhalos[0]]
	sPosy = sPosy[idofbigenoughhalos[0]]
	sPosz = sPosz[idofbigenoughhalos[0]]
	sMass = sMass[idofbigenoughhalos[0]]
	
	print(len(coslist))
	
	xlistplot = []
	ylistplot = []
	coslistplot =[]
	tol = .05e0
	
	for j in range(len(mlist)):
		xlistplot.append([])
		ylistplot.append([])
		coslistplot.append([])
		for i in range(len(HostMassRatio)):
			if fabs(HostMassRatio[i] - mlist[j]) < tol: # HostMassRatio = host/partner
				if tag[i] == 0: #A
					xlistplot[j].append((-sPosx[i] + dsep[i])/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])
				else: #B
					xlistplot[j].append(sPosx[i]/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])					
				coslistplot[j].append(coslist[i]) #only minor host distn is plotted in histogram
		
			if fabs(1e0/HostMassRatio[i] - mlist[j]) < tol:
				if tag[i] == 0:
					xlistplot[j].append((-sPosx[i] + dsep[i])/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])
				else:
					xlistplot[j].append(sPosx[i]/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])

	binno = 40
	plotpath = os.path.join(os.getcwd(),"roche_4panels.png")	

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')	
	fig = plt.figure()
	fig.set_size_inches(12.5, 2.5, forward=True)

	gsmom = gridspec.GridSpec(1, 1, hspace = .2)		
#	gsmom = gridspec.GridSpec(1, 2, hspace = 0,  height_ratios = [.6, 1]) 
	#make nested gridspecs
	gs1 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec = gsmom[0], wspace = .0, hspace = .0)
#	gs2 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec = gsmom[1], wspace = .0, hspace = .0)
	
	ax_histlist = []
	for i in range(4):
		ax_histlist.append(plt.subplot(gs1[i]))	
		
	xmin, xmax, ymin,ymax = -.5e0, 1.5e0, -.6e0, .6e0
	xlims, ylims = [xmin,xmax], [ymin,ymax]
	nxbins, nybins, nbins = binno, binno, binno
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	
	cosbins = np.arange(-1.0, 1.0 + (1.0-(-1.0))/nbins, (1.0-(-1.0))/nbins)
	coslocs = np.arange(-1.0 + .5 * (1.0-(-1.0))/nbins, 1.0 + .5*(1.0-(-1.0))/nbins, (1.0-(-1.0))/nbins)
	for i in range(len(mlist)):
		ax_histlist[i].set_xlim(-1, 1)
		ax_histlist[i].set_ylim(0.38, .88)
		ticklabels = ax_histlist[i].get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(10)
			
		#Plot the histograms
		histx, bin_edgesx = np.histogram(coslistplot[i], bins=cosbins)
		Qminus10 = sum(histx[:int(binno*0.1)])
		Qplus90 = sum(histx[int(binno*0.9):])
		
		meanAB, sigmaAB = markovN(len(coslistplot[i]))
		
		ax_histlist[i].plot(coslocs, np.array(histx)/sum(histx) * float(binno/2), color='blue', linewidth = 1.5)
		ax_histlist[i].text(-0.3, .75, "$q = %.2f, %.0f \sigma$ " % (Qplus90/Qminus10, (Qplus90/Qminus10 - meanAB)/sigmaAB), color='red', fontsize=12, fontweight='bold')	
		ax_histlist[i].set_title(" %.1f $ < M_{2}/M_{1} <$ %.1f" % (mlist[i]-0.05, mlist[i]+0.05))
		ax_histlist[i].set_xlabel("$\cos \\theta$",fontsize=14)
		
		#Make the tickmarks pretty
		if i != 0:
			ax_histlist[i].yaxis.set_ticklabels([])	

		ticklabels = ax_histlist[0].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(10)
				
		# Remove the inner axes numbers of the histograms
		ax_histlist[i].yaxis.set_major_locator(MaxNLocator(5))	
##		if i != len(mlist)-1:	
		ax_histlist[i].xaxis.set_major_locator(MaxNLocator(5, prune='upper'))
		
	ax_histlist[0].set_ylabel("$p(\cos \\theta)$",fontsize=14)
	plt.savefig(plotpath, transparent=True, dpi=500,bbox_inches='tight')
	plt.close()
	
def drawdata():
	hostcatID, coslist, HostMassRatio, covera, bovera, coverb, sToHmass, sToPmass, sMass, sPosx, sPosy, sPosz, dsep, distsat, relDistSat, tag = loadtxt(hrcosine_distn_file_roche, skiprows=1, unpack=True)
	
	toosmallhalo = (sMass /particlemass < 20)		
	idofbigenoughhalos = np.where(toosmallhalo == 0)
	
	coslist = coslist[idofbigenoughhalos[0]]
	relDistSat = relDistSat[idofbigenoughhalos[0]]
	dsep = dsep[idofbigenoughhalos[0]]
	HostMassRatio = HostMassRatio[idofbigenoughhalos[0]]
	sPosx = sPosx[idofbigenoughhalos[0]]
	sPosy = sPosy[idofbigenoughhalos[0]]
	sPosz = sPosz[idofbigenoughhalos[0]]
	sMass = sMass[idofbigenoughhalos[0]]
	
#	print(np.array(sMass[:100],dtype = int)/2.6e6)
#	print(sMass[:100])
	print(len(coslist))
	
	xlistplot = []
	ylistplot = []
	coslistplot =[]
	tol = .05e0
	
	for j in range(len(mlist)):
		xlistplot.append([])
		ylistplot.append([])
		coslistplot.append([])
		for i in range(len(HostMassRatio)):
			if fabs(HostMassRatio[i] - mlist[j]) < tol: # HostMassRatio = host/partner
				if tag[i] == 0: #A
					xlistplot[j].append((-sPosx[i] + dsep[i])/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])
				else: #B
					xlistplot[j].append(sPosx[i]/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])					
				coslistplot[j].append(coslist[i]) #only minor host distn is plotted in histogram
		
			if fabs(1e0/HostMassRatio[i] - mlist[j]) < tol:
				if tag[i] == 0:
					xlistplot[j].append((-sPosx[i] + dsep[i])/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])
				else:
					xlistplot[j].append(sPosx[i]/dsep[i])
					ylistplot[j].append(sPosy[i]/dsep[i])
#		print(len(coslistplot[j]))		

	binno = 40
	plotpath = os.path.join(os.getcwd(),"roche_4panels.png")	

	
	fig = plt.figure()
	fig.set_size_inches(10.5, 4.1, forward=True)
	mpl.rcParams['text.usetex'] = True 
	mpl.rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
	mpl.rcParams['font.family'] = 'sans-serif'
	mpl.rcParams['font.sans-serif'] = 'cm'	
		
	gsmom = gridspec.GridSpec(2, 1, hspace = .2,  height_ratios = [.6, 1]) 
	#make nested gridspecs
	gs1 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec = gsmom[0], wspace = .0, hspace = .0)
	gs2 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec = gsmom[1], wspace = .0, hspace = .0)
	
	ax_histlist = []
	for i in range(4):
		ax_histlist.append(plt.subplot(gs1[i]))	
		
	ax_list = []
	for i in range(4):
		ax_list.append(plt.subplot(gs2[i]))
							
	xmin, xmax, ymin,ymax = -.5e0, 1.5e0, -.6e0, .6e0
	xlims, ylims = [xmin,xmax], [ymin,ymax]
	nxbins, nybins, nbins = binno, binno, binno
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	
	cosbins = np.arange(-1.0, 1.0 + (1.0-(-1.0))/nbins, (1.0-(-1.0))/nbins)
	coslocs = np.arange(-1.0 + .5 * (1.0-(-1.0))/nbins, 1.0 + .5*(1.0-(-1.0))/nbins, (1.0-(-1.0))/nbins)
	for i in range(len(mlist)):
		ax_histlist[i].set_xlim(-1, 1)
		ax_histlist[i].set_ylim(0.35, 1.05)
		ticklabels = ax_histlist[i].get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(10)
			
		#Plot the histograms
		histx, bin_edgesx = np.histogram(coslistplot[i], bins=cosbins, normed = True)
		Qminus10 = sum(histx[:int(binno*0.1)])
		Qplus90 = sum(histx[int(binno*0.9):])
		
		meanAB, sigmaAB = markovN(len(coslistplot[i]))
		
		ax_histlist[i].plot(coslocs, histx, color='blue', linewidth = 1.5)
		ax_histlist[i].text(-0.5, .9, "$q = %.2f, %.0f \sigma$ " % (Qplus90/Qminus10, (Qplus90/Qminus10 - meanAB)/sigmaAB), color='red', fontsize=12, fontweight='bold')	
		ax_histlist[i].set_title(" %.1f $ < M_{2}/M_{1} <$ %.1f" % (mlist[i]-0.05, mlist[i]+0.05))
		ax_histlist[i].set_xlabel("$\cos \\theta$",fontsize=14)
		
		#Make the tickmarks pretty
		if i != 0:
			ax_histlist[i].yaxis.set_ticklabels([])	

		ticklabels = ax_histlist[0].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(10)
				
		# Remove the inner axes numbers of the histograms
		ax_histlist[i].yaxis.set_major_locator(MaxNLocator(4, prune='lower'))	
		if i != len(mlist)-1:	
			ax_histlist[i].xaxis.set_major_locator(MaxNLocator(5, prune='upper'))
		
		H, xedges, yedges = np.histogram2d(ylistplot[i], xlistplot[i], bins=(ybins, xbins))
		im = ax_list[i].imshow(log(H), extent=[xmin,xmax,ymin,ymax], cmap = 'jet', interpolation='nearest', aspect = 'equal', vmin=0, vmax=log(700))
		levels = np.arange(1, 6, 0.45)
		if i != len(mlist) - 1:
			ax_list[i].xaxis.set_major_locator(MaxNLocator(5, prune='upper'))		
		
		extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
#		ax_list[i].set_aspect(1)
		ax_list[i].contour(log(H), levels, colors=['black'], linewidths = 0.8, alpha = 0.37, extent=extent)
#		ax_list[i].xaxis.set_major_locator(MaxNLocator(5))
		ax_list[i].xaxis.set_minor_locator(AutoMinorLocator(10))
		ax_list[i].yaxis.set_major_locator(MaxNLocator(4))
#		ax_list[i].set_xlabel("$x$",fontsize=14)
		ax_list[0].set_ylabel("$y$",fontsize=14, rotation = 0)
		
		if i != 0:
			ax_list[i].yaxis.set_ticklabels([])	
		ax_list[i].set_xlabel("$x$",fontsize=14, rotation = 0)		
		ax_list[i].set_xlim(xmin, xmax)
		ax_list[i].set_ylim(ymin, ymax)
		
		ticklabels = ax_list[i].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(10)
		ticklabels = ax_list[i].get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(10)
	t = np.arange(0,8,1)
	
	plt.colorbar(im, ticks=t, format='$%.1f$',cax = fig.add_axes([0.91, 0.2, 0.02, 0.6]))
	
	ax_histlist[0].set_ylabel("$p(\cos \\theta)$",fontsize=14)

#	im.set_clim(vmin=0, vmax=10)
	plt.savefig(plotpath, transparent=True, dpi=500,bbox_inches='tight')
	plt.close()
	
def plotNvsmassratio():
	nbins = 40
	hostcatID, coslist, HostMassRatio, covera, bovera, coverb, sToHmass, sToPmass, sMass, sPosx, sPosy, sPosz, dsep, distsat, relDistSat, tag = loadtxt(hrcosine_distn_file_roche, skiprows=1, unpack=True)
	
	toosmallhalo = (sMass /particlemass < 20)		
	idofbigenoughhalos = np.where(toosmallhalo == 0)
	
	HostMassRatio = HostMassRatio[idofbigenoughhalos[0]]
	
	MRmax = 1.0
	Massratiobins = np.arange(0, MRmax + MRmax/nbins, MRmax/nbins)
		
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
	fig.savefig('N_vs_Massratio.jpg', dpi=500, bbox_inches='tight')
	
def plotlocalprince(pairno, mA, mB, sathalos, satofAm, satofBm, princeaxesAB, vAB1, hostApos, hostBpos, plotname, dsep, lim):
	plt.clf()
	fig = plt.figure(figsize=(8,8))
	ax = fig.add_subplot(221, projection='3d')

	[satofA, satofB] = sathalos
	m1, m2, m3 = middlecoordmpc(hostApos, hostBpos)
	ax.set_xlim3d(m1 - dsep * lim, m1 + dsep * lim)
	ax.set_ylim3d(m2 - dsep * lim, m2 + dsep * lim)
	ax.set_zlim3d(m3 - dsep * lim, m3 + dsep * lim)
		
	for i in range(len(satofA)):
		mass = satofAm[i]
		ax.scatter(satofA[i][0], satofA[i][1], satofA[i][2], zdir='z', color = 'red', s = mass/plotmassscale)

	for i in range(len(satofB)):
		mass = satofBm[i]
		ax.scatter(satofB[i][0], satofB[i][1], satofB[i][2], zdir='z', color = 'blue', s = mass/plotmassscale)

	shiftedaxesA = [[],[],[]]
	shiftedaxesB = [[],[],[]]
	for i in range(3):
		for j in range(3):
			shiftedaxesA[i].append(princeaxesAB[i][j]/5e0) #5 for scaling length of the eigenvector
			shiftedaxesB[i].append(princeaxesAB[i][j]/5e0)

	for i in range(3):
#		ax.plot([dsep, shiftedaxesB[i][0]+dsep], [0, shiftedaxesB[i][1]], [0, shiftedaxesB[i][2]], color = 'green')
		ax.plot([hostApos[0], hostApos[0] + shiftedaxesB[i][0]], [hostApos[1], hostApos[1] + shiftedaxesB[i][1]], [hostApos[2], hostApos[2] + shiftedaxesB[i][2]], color = 'green')
	
	ax.plot([hostApos[0], hostApos[0]+vAB1[0]], [hostApos[1], hostApos[1]+vAB1[1]], [hostApos[2], hostApos[2]+vAB1[2]], color = 'purple')
	
	#take care of periodic boundary condition for binary
#	hostApos, hostBpos, marker = coordwrap(hostApos, hostBpos) #if marker[i] = 1, hostApos 63.5 -> -.5; marker[i] = 2 hostBpos 63.5 -> -.5
	ax.scatter(hostApos[0], hostApos[1], hostApos[2], alpha = 0.3, color = 'black', s = mA/plotmassscale/10) 	#halo size proportional to mass
	ax.scatter(hostBpos[0], hostBpos[1], hostBpos[2], alpha = 0.3, color = 'blue', s = mB/plotmassscale/10)
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("z")
	
	ax2 = fig.add_subplot(222)
	ax2.set_xlim(m1-dsep*lim, m1+dsep*lim)
	ax2.set_ylim(m2-dsep*lim, m2+dsep*lim)
		
	for i in range(len(satofA)):
		mass = satofAm[i]
		ax2.scatter(satofA[i][0], satofA[i][1], color = 'red', s = mass/plotmassscale)

	for i in range(len(satofB)):
		mass = satofBm[i]
		ax2.scatter(satofB[i][0], satofB[i][1], color = 'blue', s = mass/plotmassscale)
		
	ax2.scatter(hostApos[0], hostApos[1], color = 'red', s = mA/plotmassscale/10) 	#halo size proportional to mass
	ax2.scatter(hostBpos[0], hostBpos[1], color = 'blue', s = mB/plotmassscale/10)
	ax2.plot([hostApos[0], hostApos[0]+vAB1[0]], [hostApos[1], hostApos[1]+vAB1[1]], color = 'purple')
	ax2.set_xlabel("x")
	ax2.set_ylabel("y")
	
	colorlist = ['purple', 'green', 'blue']

	for i in range(3):
#		ax2.plot([0, shiftedaxesA[i][0]], [0, shiftedaxesA[i][1]], color = colorlist[i])
#		ax2.plot([dsep, shiftedaxesB[i][0]+dsep], [0, shiftedaxesB[i][1]], color = colorlist[i])
		ax2.plot([hostApos[0], hostApos[0]+shiftedaxesB[i][0]], [hostApos[1],hostApos[1]+shiftedaxesB[i][1]], color = 'green')
		
	ax3 = fig.add_subplot(223)
	ax3.set_xlim(m2-dsep*lim, m2+dsep*lim)	
	ax3.set_ylim(m3-dsep*lim, m3+dsep*lim)
	ax3.set_xlabel("y")
	ax3.set_ylabel("z")	
	for i in range(len(satofA)):
		mass = satofAm[i]
		ax3.scatter(satofA[i][1], satofA[i][2], color = 'red', s = mass/plotmassscale)

	for i in range(len(satofB)):
		mass = satofBm[i]
		ax3.scatter(satofB[i][1], satofB[i][2], color = 'blue', s = mass/plotmassscale)
	ax3.scatter(hostApos[1], hostApos[2], color = 'red', s = mA/plotmassscale/10) 	#halo size proportional to mass
	ax3.scatter(hostBpos[1], hostBpos[2], color = 'blue', s = mB/plotmassscale/10)
	ax3.plot([hostApos[1], hostApos[1]+vAB1[1]], [hostApos[2], hostApos[2]+vAB1[2]], color = 'purple')	

	for i in range(3):
#		ax3.plot([0, shiftedaxesA[i][1]], [0, shiftedaxesA[i][2]], color = colorlist[i])
#		ax3.plot([0, shiftedaxesB[i][1]], [0, shiftedaxesB[i][2]], color = colorlist[i])
		ax3.plot([hostApos[1], hostApos[1]+shiftedaxesB[i][1]], [hostApos[2],hostApos[2]+shiftedaxesB[i][2]], color = 'green')	

	ax4 = fig.add_subplot(224)
	ax4.set_xlim(m1-dsep*lim, m1+dsep*lim)
	ax4.set_ylim(m3-dsep*lim, m3+dsep*lim)

	ax4.set_xlabel("x")
	ax4.set_ylabel("z")
		
	for i in range(len(satofA)):
		mass = satofAm[i]
		ax4.scatter(satofA[i][0], satofA[i][2], color = 'red', s = mass/plotmassscale)

	for i in range(len(satofB)):
		mass = satofBm[i]
		ax4.scatter(satofB[i][0], satofB[i][2], color = 'blue', s = mass/plotmassscale)
	ax4.scatter(hostApos[0], hostApos[2], color = 'red', s = mA/plotmassscale/10) 	#halo size proportional to mass
	ax4.scatter(hostBpos[0], hostBpos[2], color = 'blue', s = mB/plotmassscale/10)		
	ax4.plot([hostApos[0], hostApos[0]+vAB1[0]], [hostApos[2], hostApos[2]+vAB1[2]], color = 'purple')	
	for i in range(3):
		ax4.plot([hostApos[0], hostApos[0]+shiftedaxesB[i][0]], [hostApos[2],hostApos[2]+shiftedaxesB[i][2]], color = 'green')
#		ax4.plot([0, shiftedaxesA[i][0]], [0, shiftedaxesA[i][2]], color = colorlist[i])
#		ax4.plot([dsep, shiftedaxesB[i][0]+dsep], [0, shiftedaxesB[i][2]], color = colorlist[i])	
	fig.savefig(plotname, dpi=500)
	plt.close()
	
				
if __name__ == "__main__":

#	Rockstarfilelist0 = []
#	Rockstar_files = os.listdir(Rockstar_hlist_address)
#	
#	for i in range(len(Rockstar_files)):
#		filename = Rockstar_files[i]
#		if filename.endswith("list"):
#			Rockstarfilelist0.append(filename)
#	Rockstarfilelist0.sort(key = natural_keys)
#	filename_high_res = Rockstarfilelist0[::-1][0]
#		
#	# copy hlist file to home address
#	snapnum = filename_high_res[6:12]
#	datfoldr = "snapshot" + snapnum
#		
#	filepath_catalog = Rockstar_hlist_address + filename_high_res
#	newfilepath = datfoldr + "/" + filename_high_res
#	os.chdir("./" + datfoldr)
#	
	mlist = [.05e0, .15e0, .25e0, .35e0]
#	drawRochePot()
#	sys.exit()
#	hrLGcathaloID, hrLGm, hrLGx, hrLGy, hrLGz, hrLGRvir, hrLGb, hrLGc = col_reader5(high_res_trueLG_pop, IDcol, masscol, xcol, ycol, zcol, Rvircol, bcol, ccol)
#	cosdistn(hrLGcathaloID, hrLGm, hrLGx, hrLGy, hrLGz, hrLGRvir, hrLGb, hrLGc)

#	drawdata()
#	drawdata_onepanel()
#	drawdata_onerow()

