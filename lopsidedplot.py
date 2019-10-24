import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from constant import *
from helpers import *
from merger_sat_analysisplot import *
from density2dplot import *

def plot_cosine_hist_wsigma(binno):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
	fig = plt.figure()
	fig.set_size_inches(10, 13, forward=True)
	ax = fig.add_subplot(211)
	
	coslist = loadtxt(hrcosine_distn_file, usecols = (0,), skiprows=1, unpack=True)

	idofbigenoughhalos = np.loadtxt('idofbigenoughhalos_z=0.dat')
	idofbigenoughhalos = np.array(idofbigenoughhalos, dtype = int)

	coslist = coslist[idofbigenoughhalos]
	
	indcospos = np.where(np.array(coslist) > 0)
	indcosneg = np.where(np.array(coslist) < 0)	
	
#	print(mean(coslist[indcospos[0]]), -mean(coslist[indcosneg[0]]))
#	print(markov(len(coslist)))
#	sys.exit()
	
	indpos = np.where(np.array(coslist) > 0.8)
	indneg = np.where(np.array(coslist) < -0.8)	
	
	Qplus90 = len(indpos[0])
	Qminus10 = len(indneg[0])

	meanAB, sigmaAB = markovN(len(coslist))
	print(meanAB, sigmaAB)
	
	meanabs, sigmaabs = markovNabs(len(coslist))
	print(meanabs, sigmaabs)
	
	hist, bin_edges = np.histogram(coslist, binno * 2)
	xlist = np.arange(-1e0 + 2.5e-2, 1e0, 5e-2)

#	for tl in ax.get_xticklabels():
#		tl.set_color('b')
    		
	ticklabels = ax.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
	ticklabels = ax.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)	
	
	hlist = np.array(hist,dtype = float)/float(sum(hist)) * float(binno)
	ax.plot(xlist, hlist, linewidth=2, color='b')
	ax.set_xlabel('$\cos\\theta$',fontsize=20)
	ax.set_ylabel('$p(\\cos\\theta)$',fontsize=20)
	
#	ax.spines['bottom'].set_color('blue') 
	
	cosmean = np.median(coslist)
#	cossigma = np.std(coslist)
	ax.axvline(cosmean, linewidth=1, color='b', alpha = 0.7) #mean of cosine	
#	print (cosmean, cossigma) #0.00175693856252 0.602460001536

	ax.text(-.5, 0.66, "$q=%.3f, %.0f \sigma$" % (Qplus90/Qminus10, abs((Qplus90/Qminus10 - meanAB)/sigmaAB)), color='blue', fontsize=20, fontweight='bold')
#	legend = ax.legend(loc="upper left", markerscale=4., shadow=True, fontsize=16)
	ax.set_ylim(4.2e-1, 6.8e-1)
	ax.set_xlim(-1e0, 1e0)

	poscoslist = []
	negcoslist = []
	for i in coslist:
		if i > 0:
			poscoslist.append(i)
		else:
			negcoslist.append(abs(i))
				
	plotREFname = os.path.join(os.getcwd(), "paper", "coshist_wsigma20_reflected.pdf")
		
	ax2 = fig.add_subplot(212)
	ax2.set_ylim(4.2e-1, 6.8e-1)
	ax2.set_xlabel('$|\cos\\theta|$',fontsize=20)
	ax2.set_ylabel('$p(\\cos\\theta)$',fontsize=20)
		
	ax2.plot(xlist[binno:2*binno], hlist[binno:2*binno], 'b-', linewidth=2, label = "$\\cos \\theta > 0$")
	ax2.fill_between(xlist[binno:2*binno], hlist[binno:2*binno], hlist[:binno][::-1], hatch = '\\', facecolor="none", edgecolor="grey", linewidth=0.0)
	ax2.plot(xlist[:binno][::-1]+1, hlist[:binno], 'b--', linewidth=2, label = "$\\cos \\theta < 0$")

	ticklabels = ax2.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(18)
	ticklabels = ax2.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(18)	
			
	poscosmean = np.median(poscoslist)
	negcosmean = np.median(negcoslist)	

	ax2.axvline(poscosmean, linewidth=1, color='b',linestyle = '-', alpha = 0.7) #mean of cosine	
	ax2.axvline(negcosmean, linewidth=1, color='b',linestyle = '--', alpha = 0.7) #mean of cosine
	
	ax2.set_xlim(0, 1e0)
#	ax2.spines['top'].set_color('red') 
	legend = ax2.legend(loc=(0.02,0.76), markerscale=4., shadow=True,fontsize=18)

#	for tl in ax2.get_xticklabels():
#		tl.set_color('r')
			
	fig.savefig(plotREFname, dpi=500, bbox_inches='tight')

def densitylinesplot(binno, xlabel, ylabel, plotname):
	coslist, distsat = loadtxt(hrcosine_distn_file, usecols = (0, 10), skiprows=1, unpack=True)	
	print(len(coslist))
	
	idofbigenoughhalos = np.loadtxt('idofbigenoughhalos_z=0.dat')
	idofbigenoughhalos = np.array(idofbigenoughhalos, dtype = int)
	
	cos = coslist[idofbigenoughhalos]
	r = distsat[idofbigenoughhalos]
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	dlims = [0.1,0.2,0.3,0.4,0.5]
	coslists = []
	for j in range(len(dlims)-1):
		coslists.append([])
	for i in range(len(r)):
		for j in range(len(dlims)-1):
			if dlims[j] < r[i] < dlims[j+1]:
				coslists[j].append(cos[i])
	
	plt.clf()
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xlabel(xlabel,fontsize=15)
	ax.set_ylabel(ylabel,fontsize=15)
	fig.set_size_inches(5, 4.0, forward=True)
	xlist = np.arange(-1e0 + 5e-2, 1e0, 1e-1)

	linestyles = ['-', '--', '-.', ':']
	
	for j in range(len(coslists))[::-1]:
		coslist = coslists[j]

		indpos = np.where(np.array(coslist) > 0.8)
		indneg = np.where(np.array(coslist) < -0.8)	
		
		Qplus90 = len(indpos[0])
		Qminus10 = len(indneg[0])
		
		meanAB, sigmaAB = markovN(len(coslist))
		print(meanAB, sigmaAB)
		
		hist, bin_edges = np.histogram(coslist, binno)
		
		q = float(Qplus90)/float(Qminus10)
		hlist = np.array(hist,dtype = float)/float(sum(hist)) * float(binno) /2e0
		ax.plot(xlist, hlist, linestyle=linestyles[j], label = " $%.1f < d_{\mathrm{sat}}/d_{\mathrm{sep}} < %.1f, q = %.3f, %.0f \sigma$ " % (dlims[j], dlims[j + 1], q, (q - meanAB) / sigmaAB) )
		
		cosmean = np.mean(coslist)
		ax.set_ylim(0.4,0.82)
		ax.set_xlim(-1,1)
		
	legend = ax.legend(loc="upper center", shadow=True, fontsize = 12)
	fig.savefig(plotname, dpi=500, bbox_inches='tight')
		
def rtheta_densityplot():
	#plot distance vs. angle for small satellites when a feature is low or high
	coslist, HostMassRatio, covera, bovera, coverb, sToHmass, sToPmass, sMass, dsep, distsat, relDistSat = loadtxt(hrcosine_distn_file, skiprows=1, unpack=True)

	idofbigenoughhalos = np.loadtxt('idofbigenoughhalos_z=0.dat')
	idofbigenoughhalos = np.array(idofbigenoughhalos, dtype = int)

	HostMassRatio = HostMassRatio[idofbigenoughhalos]
	coslist = coslist[idofbigenoughhalos]
	relDistSat = relDistSat[idofbigenoughhalos]

	filename = ["hostRelMass", "cOvera", "satMass", "dist_of_bin"]
	featname = [r"$\boldsymbol{\frac{M_{1}}{M_{2}}}$", "c/a", r"$\boldsymbol{M_{sat}}$", r"$\boldsymbol{d_{\mathrm{sep}}}$"]
	
	index = np.where(HostMassRatio < 1.)
	
	featlist = [HostMassRatio[index[0]], covera[index[0]], sMass[index[0]], dsep[index[0]]]
	
	coslist = coslist[index[0]]
	relDistSat = relDistSat[index[0]]
	
	#plot r vs costheta for data biased on high or low features
	xlistplot = []
	ylistplot = []		
	for j in range(len(featlist)):
		sort_tup_list = []
		for i in range(len(coslist)):
			sort_tup_list.append((coslist[i], featlist[j][i], relDistSat[i]))

		sort_tup_list.sort(key=lambda tup: tup[1]) #sort according to jth feature

		sorted_distlist = []
		for i in range(len(sort_tup_list)):
			sorted_distlist.append(sort_tup_list[i][2])

		sorted_coslist = []
		for i in range(len(sort_tup_list)):
			sorted_coslist.append(sort_tup_list[i][0])
			
		n1 = len(sorted_distlist)
	
		chunk_num = 10e0 #10% chunks
		chunk_size = int(n1 / chunk_num)
		groups = []
		
		for i in range(0, n1, chunk_size):
			groups.append([sorted_coslist[i : i + chunk_size], sorted_distlist[i : i + chunk_size]])

		group_tile = [0,9] #top and bottom 10%
		xlistplot.append([])
		ylistplot.append([])
		for i in group_tile:
			xlistplot[j].append(groups[i][0])
			ylistplot[j].append(groups[i][1])
			
	plotpath = os.path.join(os.getcwd(), "paper", "density_r_vs_cos_8panels_minorhost.pdf")				
	densityplot_8panels(xlistplot, ylistplot, 'cos $\\theta$', '$d_{\mathrm{sat}}/d_{\mathrm{sep}}$', featname, plotpath, [-1e0,1e0], 40)			

def rtheta_densityplot_all():
	#plot distance vs. angle for small satellites when a feature is low or high
	coslist, HostMassRatio, covera, bovera, coverb, dsep, sMass, satx, saty, satz, relDistSat = loadtxt(hrcosine_distn_file, skiprows=1, unpack=True)

#	[cosangle, mass of satellite, dist to host, satx, saty, satz]
#	f.write(str(cosA[i][0]) + "\t" + str(mA/mB) + "\t" + str(cA) + "\t" + str(bA) + "\t" + str(cA/bA) + "\t" + str(Sepdis) + "\t" + str(cosA[i][1]) + "\t" + str(cosA[i][3]) + "\t" + str(cosA[i][4]) + "\t" + str(cosA[i][5]) + "\t" + str(cosA[i][2]/Sepdis) + "\n")
	
	idofbigenoughhalos = np.loadtxt('idofbigenoughhalos_z=0.dat')
	idofbigenoughhalos = np.array(idofbigenoughhalos, dtype = int)

	HostMassRatio = HostMassRatio[idofbigenoughhalos]
	coslist = coslist[idofbigenoughhalos]
	relDistSat = relDistSat[idofbigenoughhalos]
	
#	print(len(coslist))
#	sys.exit()
	index1 = np.where(HostMassRatio < 1.)
	index2 = np.where(HostMassRatio > 1.)
	
	coslist1 = coslist[index1[0]]
	relDistSat1 = relDistSat[index1[0]]

	coslist2 = coslist[index2[0]]
	relDistSat2 = relDistSat[index2[0]]
		
	plotpath = os.path.join(os.getcwd(), "paper", "density_r_vs_cos_minorandmajor_host.pdf")				
	densityplot_2panels(coslist, relDistSat, coslist2, relDistSat2, '$\cos \\theta$', '$d_{\mathrm{sat}}/d_{\mathrm{sep}}$', plotpath, [-1e0,1e0], [0,.5e0], [0.4,0.75], 2.5, 40)

def MassratioHostCumulative(cumtag):
	TLGID, Tlgbm, Tlgbx, Tlgby, Tlgbz, TlgRvir, Tlgb, Tlgc = col_reader5(hrtrueLG_wosubhalo_biggap_filtered, IDcol, masscol, xcol, ycol, zcol, Rvircol, bcol, ccol)
	print(TlgRvir[:100])
	largeRvir = np.where(TlgRvir> 100)
	print(len(largeRvir[0]), len(TlgRvir))
	sys.exit()
	MRlist = []
	for i in range(len(Tlgbm)):
		if i % 2==0:
			MR = Tlgbm[i]/Tlgbm[i+1]
#			print(Tlgbm[i],Tlgbm[i+1])
			if MR > 1:
				MR = 1e0/MR
			MRlist.append(MR) 
	print(len(MRlist))
	MRmax = 1e0
	nbins = 40
	Massratiobins = np.arange(0, MRmax + MRmax/nbins, MRmax/nbins)
	
	histx, bin_edgesx = np.histogram(MRlist, bins=Massratiobins)
#	cumulative = np.cumsum(histx)

	cumulative = []
	for i in range(len(histx)):
		cumulative.append(sum(histx[:i+1]))
	
	MRlocs = np.arange(0 + .5 * MRmax/nbins, MRmax + .5*MRmax/nbins, MRmax/nbins)
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')	
	fig = plt.figure()
	fig.set_size_inches(8, 5, forward=True)
	ax = fig.add_subplot(111)
	if cumtag ==1:	
		ax.plot(MRlocs, cumulative/sum(histx), color='blue', linewidth = 1.5)
	else:
		ax.plot(MRlocs, histx/sum(histx), color='blue', linewidth = 1.5)
		
	ax.set_xlim(0, MRmax)	
	ax.set_xlabel('$M_2/M_1$', fontsize=20)
	if cumtag ==1 :
		ax.set_ylabel('$N_{pair}(<M_2/M_1)$',fontsize=20)
	else:
		ax.set_ylabel('$N_{pair}$',fontsize=20)
	plt.grid(True)
	if cumtag ==1:
		fig.savefig('Nhost_vs_Massratio_cumulative.jpg', dpi=500, bbox_inches='tight')
	else:
		fig.savefig('Nhost_vs_Massratio.jpg', dpi=500, bbox_inches='tight')
if __name__ == "__main__":

	sMass = loadtxt(hrcosine_distn_file, usecols = (6,), skiprows=1, unpack=True)
	print(min(sMass), min(sMass)/2.6013e6)
#	bigenoughhalo = (sMass /2.6013e6> 20)
#	print(min(sMass))		
#	idofbigenoughhalos = np.where(bigenoughhalo == 1)
#	print(len(idofbigenoughhalos[0]))
#	np.savetxt('idofbigenoughhalos_z=0.dat' , idofbigenoughhalos[0], delimiter=" ")		

#	MassratioHostCumulative(0)	
#	MassratioHostCumulative(1)		
#	sys.exit()
#	
#	plot_cosine_hist_wsigma(20)
##	
#	sys.exit()

#	densitylinesplot(20, ' $\cos\\theta$', '$p(\cos\\theta)$', os.path.join(os.getcwd(), "paper", "rtheta_rdep_density.pdf"))
	
#	sys.exit()
	
	rtheta_densityplot_all()
