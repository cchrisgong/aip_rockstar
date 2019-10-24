from numpy import loadtxt, histogram
import numpy as np
import sys
from constant import *
import os
from helpers import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
import matplotlib as mpl
from numpy import amax, amin, median, sin, cos, arccos, std, log, mean, log10, log2
from mpl_toolkits.mplot3d import Axes3D
from density2dplot import *
from merger_analysis import plotread_host

def myhist(var, pmax, pmin, binnum):
	binlist = np.arange(pmin, pmax + (pmax - pmin)/binnum, (pmax - pmin)/float(binnum))
	xlocs = np.arange(pmin + .5 * (pmax - pmin)/float(binnum), pmax + .5 * (pmax - pmin)/float(binnum), (pmax - pmin)/float(binnum))
#	print(binlist, np.where(np.isnan(var)), len(var))
	hist, bin_edgesx = np.histogram(var, bins=binlist)
	print(pmin, pmax, binnum, xlocs)
	return xlocs, hist/sum(hist) * float(binnum) /2e0

def percentile(var_list, i):
	var_list = np.sort(var_list)
	list_perc = [.68,.95, .997]
	N = len(var_list)
	upper_perc = .5 + list_perc[i-1] / 2e0
	lower_perc = .5 - list_perc[i-1] / 2e0
	return var_list[round(lower_perc * N)], var_list[round(upper_perc * N)]

def plotread(filelist0, l, variabletag):
	idofbigenoughhalosA = np.loadtxt('idofbigenoughhalosA.dat')
	idofbigenoughhalosB = np.loadtxt('idofbigenoughhalosB.dat')
#	FILE: AID_sat_1z, Acosangle_sat_1z, Av_relhost_cosangle_sat_1z, Av_relcent_cosangle_sat_1z, Apv_relhost_cosangle_sat_1z, Apv_relcent_cosangle_sat_1z, Adist_sat_1z, Am_sat_1z, Av_relhost_sat_1z

	if variabletag == 'ID':
		var_list = loadtxt(filelist0[l], usecols = (0,), unpack = 1)
		return var_list
		
	elif variabletag == 'descID':
		var_list = loadtxt(filelist0[l], usecols = (1,), unpack = 1)
		return var_list
		
	elif variabletag == 'logmass_sat':
		m_sat_list_1z = loadtxt(filelist0[l], usecols = (8,), unpack = 1)
		pmax = 12.
#		pmin = min(m_sat_list_1z)
		var_list = np.log10(m_sat_list_1z)
	
	elif variabletag == 'sat_cosangle_pos':
		var_list = loadtxt(filelist0[l], usecols = (2,), unpack = 1)
		pmax = 1e0
		pmin = -1e0
#		var_list = np.array([x for x in cosangle_sat_1z if not math.isnan(x)])

	elif variabletag == 'sat_sinangle_pos':
		var_list = loadtxt(filelist0[l], usecols = (10,), unpack = 1) 
		pmax = 1e0
		pmin = -1e0
#		var_list = np.array([x for x in cosangle_sat_1z if not math.isnan(x)])
	
	elif variabletag == 'sat_cosangle_vel_relhost':
		var_list = loadtxt(filelist0[l], usecols = (3,), unpack = 1)
		pmax = 1e0
		pmin = -1e0

	elif variabletag == 'sat_cosangle_vel_relcent':
		var_list = loadtxt(filelist0[l], usecols = (4,), unpack = 1)
		pmax = 1e0
		pmin = -1e0
				
	elif variabletag == 'dist_sat': #normalized by dsep0
		var_list = loadtxt(filelist0[l], usecols = (7,), unpack = 1)
		pmax = 0.8e0
		pmin = 0e0

	elif variabletag == 'pv_relhost_cosangle':
		var_list = loadtxt(filelist0[l], usecols = (5,), unpack = 1)
		pmax = 1e0
		pmin = -1e0

	elif variabletag == 'pv_relcent_cosangle':
		var_list = loadtxt(filelist0[l], usecols = (6,), unpack = 1)
		pmax = 1e0
		pmin = -1e0
		
	elif variabletag == 'vsat':
		z = z_redshift_list[l]
		a = 1. / (z + 1.)
		var_list  = loadtxt(filelist0[l], usecols = (9,), unpack = 1)
		var_list = var_list * sqrt(a)
		pmax = 400
		pmin = 0
	
#	return pmax, pmin, var_list
	
	if filelist0[0][-5] == 'A':
		if variabletag == 'logmass_sat':
			return max(var_list[np.array(idofbigenoughhalosA, dtype = int)]), min(var_list[np.array(idofbigenoughhalosA, dtype = int)]), var_list[np.array(idofbigenoughhalosA, dtype = int)]
		else:
			return pmax, pmin, var_list[np.array(idofbigenoughhalosA, dtype = int)]
	else:
		if variabletag == 'logmass_sat':
			return max(var_list[np.array(idofbigenoughhalosB, dtype = int)]), min(var_list[np.array(idofbigenoughhalosB, dtype = int)]), var_list[np.array(idofbigenoughhalosB, dtype = int)]
		else:
			return pmax, pmin, var_list[np.array(idofbigenoughhalosB, dtype = int)]
			
def merger_sat_analysis_plot(filelist0, variabletag, fileext):
	plot_list_mat = []
	plot_list = []
	lower_sigma_bounds = []
	upper_sigma_bounds = []
	plot_redshift_list = [z_redshift_list[0]]
	
	counter = 1
	var_list0 = np.array([])
	
	for l in range(len(filelist0)):
		if l == 0:
			pmax, pmin, var_list = plotread(filelist0,l,variabletag)
			equaltoresolutionof20particles = np.where(var_list == pmin)
#			print(len(largerthanresolutionof20particles[0]), len(var_list))
			
		if z_redshift_list[l] < z_redshift_list_linear[counter]:
			pmax, pmin, var_list = plotread(filelist0,l,variabletag)
#			print(l, 10e0**pmin, filelist0)
#			largerthanresolutionof20particles = np.where(var_list > pmin)
#			print(len(largerthanresolutionof20particles[0]), len(var_list))

			print(filelist0[l], 10e0 ** var_list[equaltoresolutionof20particles[0][:20]])
#			sys.exit()
			var_list0 = np.append(var_list0, var_list)
		else:
			var_list0.flatten()
			xlocs, histx = myhist(var_list0, pmax, pmin, 100.)
			
			plot_list_mat.append(histx[::-1])
			plot_list.append([amax(var_list), amin(var_list), median(var_list), std(var_list)])	
			lower_1sigma, upper_1sigma = percentile(var_list, 1)
			lower_2sigma, upper_2sigma = percentile(var_list, 2)
			lower_3sigma, upper_3sigma = percentile(var_list, 3)
			lower_sigma_bounds.append([lower_1sigma, lower_2sigma, lower_3sigma])
			upper_sigma_bounds.append([upper_1sigma, upper_2sigma, upper_3sigma])

			counter += 1
			var_list0 = np.array([])
			if l != len(filelist0)-1:
				plot_redshift_list.append(z_redshift_list[l])
				pmax, pmin, var_list = plotread(filelist0,l,variabletag)
				var_list0 = np.append(var_list0, var_list)
	lower_sigma_bounds = np.transpose(lower_sigma_bounds)
	upper_sigma_bounds = np.transpose(upper_sigma_bounds)
	xmin = 0
	xmax = z_redshift_list[-1]
	
	fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(12,6))
	ax1.set_xlabel("z",fontsize=18)
	ax1.set_xlim(xmin, xmax)
	ax1.set_ylabel(ylabel[varname_list.index(variabletag)], fontsize=18)	
	ax1.grid(True, linestyle= '-', which= 'major',color= '0.75')
	ax1.grid(True, linestyle= '-', which= 'minor',color= '0.75')
	ax1.grid(True, which='both')
	ax1.minorticks_on()	
	plot_list = np.transpose(plot_list)
	
	ax1.plot(plot_redshift_list, plot_list[0])
	ax1.plot(plot_redshift_list, plot_list[1])
	ax1.fill_between(plot_redshift_list, plot_list[1], plot_list[0], facecolor='blue', alpha=0.2) # between max and min
	ax1.plot(plot_redshift_list, plot_list[2])
	ax1.fill_between(plot_redshift_list, lower_sigma_bounds[0], upper_sigma_bounds[0], facecolor='yellow', alpha=0.5) # one sigma
	ax1.fill_between(plot_redshift_list, lower_sigma_bounds[1], upper_sigma_bounds[1], facecolor='green', alpha=0.5) # two sigma
	ax1.fill_between(plot_redshift_list, lower_sigma_bounds[2], upper_sigma_bounds[2], facecolor='red', alpha=0.5) # three sigma
	
	ax2.set_xlabel("z", fontsize=18)	
	plot_list_mat = np.transpose(plot_list_mat)
	im = plt.imshow(plot_list_mat, extent=[xmin, xmax, pmin, pmax], cmap=plt.cm.BuPu_r, aspect = 'auto', interpolation='nearest')
	t = np.arange(amin(plot_list_mat.flatten()), amax(plot_list_mat.flatten()), (amax(plot_list_mat.flatten()) - amin(plot_list_mat.flatten())) / 5e0)
	im.set_clim(vmin=amin(plot_list_mat.flatten()), vmax=amax(plot_list_mat.flatten()))
	fig.colorbar(im, ticks=t, format='$%.2f$',cax = fig.add_axes([0.91, 0.15, 0.02, 0.7]))
	fig.savefig(os.path.join(os.getcwd(), "paper", 'traceMaxProg_' + variabletag +'_' + fileext + '.pdf'), dpi=500)
	plt.close()

def norm(val):
    return (val - z_redshift_list[0]) / float(z_redshift_list[-1] - z_redshift_list[0])
    	
def merger_sat_analysis_plot3d(variabletag, sampletag):
	pmax = 1e0
	pmin = -1e0
	
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
		
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ymin = 0
	ymax = z_redshift_list[-1]
	
	z_plotlist = [0.0, 0.5179, 0.9489, 1.5025, 1.9078,2.5511]
	
	for l in range(len(z_redshift_list)):
		if z_redshift_list[l] in z_plotlist: # over z
			
			if sampletag == 'A':
				pmax, pmin, var_list = plotread(filelist0A, l, variabletag)
			elif sampletag == 'B':
				pmax, pmin, var_list = plotread(filelist0B, l, variabletag)
			elif sampletag == 'AB':
				pmax, pmin, var_listA = plotread(filelist0A, l, variabletag)
				pmax, pmin, var_listB = plotread(filelist0B, l, variabletag)
				var_list = np.concatenate((var_listA, var_listB))
#				print(len(np.where(var_list > 0)[0])/len(var_list))
#				sys.exit()
			
#			if variabletag == 'sat_cosangle_pos':
#				indpos = np.where(var_list > 0.8)
#				indneg = np.where(var_list < -0.8)
#				Qplus90 = len(indpos[0])
#				Qminus10 = len(indneg[0])
#				meanAB, sigmaAB = markovN(len(var_list))
#				print(len(var_list))
#				print(z_redshift_list[l], Qplus90/Qminus10, (Qplus90/Qminus10 - meanAB)/sigmaAB)
##				qlist = [1.05, 1.30, 1.42, 1.59, 1.70, 1.81]
##				ax.plot(z_plotlist, qlist, zdir='x')

			x, histx = myhist(var_list, pmax, pmin, 40)
			DAT = np.column_stack((x, histx))
			np.savetxt('z=%.4f_hist_merger_' % z_redshift_list[l] + variabletag + sampletag + '.dat' , DAT, delimiter=" ")		

	for l in range(len(z_redshift_list)):
		if z_redshift_list[l] in z_plotlist: # over z
			x, histx = loadtxt('z=%.4f_hist_merger_' % z_redshift_list[l] + variabletag + sampletag + '.dat',  unpack = 1)
			print(x)
			y = np.ones(x.size) * z_redshift_list[l]
			z = histx
			ax.plot(x, y, z, color=plt.cm.gnuplot(norm(y[0])), zorder = l)
			if variabletag == 'sat_cosangle_pos':
				ax.add_collection3d(plt.fill_between(x, 0.3, z, color=plt.cm.gnuplot(norm(y[0])), alpha=0.1), zs=z_redshift_list[l], zdir='y')
			else:
				ax.add_collection3d(plt.fill_between(x, 0, z, color=plt.cm.gnuplot(norm(y[0])), alpha=0.1), zs=z_redshift_list[l], zdir='y')
				
	if variabletag == 'sat_cosangle_pos':
		ax.set_zlim(0.3, 1.)
		ax.set_zlabel("$p(\cos\\theta)$",fontsize=18)
	else:
		ax.set_zlim(0, 2.9)
		ax.set_zlabel("$p(\cos\\theta_{\mathrm{rv}})$",fontsize=18)
	ax.set_xlim(-1, 1)
	ax.set_ylim(0, 2.56)
	
	ax.set_ylabel("$z$",fontsize=18)
	ax.set_xlabel(ylabel2[varname_list2.index(variabletag)],fontsize=18)
	ax.view_init(15,120)
	
	plt.tight_layout()
	fig.savefig(os.path.join(os.getcwd(), "paper", 'traceMaxProg_' + variabletag + '3d_stacked3'+sampletag+'.pdf'), dpi=500)
	plt.close()

def merger_sat_density_plot(filelistA, filelistB, plottag):
	plot_list = []
	var_list0 = np.array([])
	z_plotlist = [0.0, 0.5179, 1.5025, 2.5511]
#	z_plotlist = [0, 0.1245, 0.3734, 2.5511]
	
	xlistall = []
	ylistall = []
	xlims = []
	ylims = []
	xnames = []
	ynames = []
	counter = 0
	redshiftnames = []
	
	for l in range(27): # over z  
		if z_redshift_list[l] in z_plotlist: # over z
			xlistall.append([])
			ylistall.append([])
			counter +=1
			redshiftnames.append('z = ' + str(round(z_redshift_list[l],1)))

			for i in range(len(variabletag1_list)): # over variables to plot
				if plottag == 'AB':
					pmax1, pmin1, var_list1A = plotread(filelistA, l, variabletag1_list[i]) # massiveornot = 0 for massive, 1 for less massive
					pmax2, pmin2, var_list2A = plotread(filelistA, l, variabletag2_list[i])
					pmax1, pmin1, var_list1B = plotread(filelistB, l, variabletag1_list[i]) # massiveornot = 0 for massive, 1 for less massive
					pmax2, pmin2, var_list2B = plotread(filelistB, l, variabletag2_list[i])						
					var_list1 = np.concatenate((var_list1A, var_list1B))
					var_list2 = np.concatenate((var_list2A, var_list2B))
				else:
					if plottag == 'A':
						filelist = filelistA
					if plottag == 'B':
						filelist = filelistB
					pmax1, pmin1, var_list1 = plotread(filelist, l, variabletag1_list[i]) # massiveornot = 0 for massive, 1 for less massive
					pmax2, pmin2, var_list2 = plotread(filelist, l, variabletag2_list[i])
				xlistall[counter - 1].append(var_list1)
				ylistall[counter - 1].append(var_list2)
				if l == 0:
					xlims.append([pmin1, pmax1])
					ylims.append([pmin2, pmax2])
					xnames.append(ylabel[varname_list.index(variabletag1_list[i])])
					ynames.append(ylabel[varname_list.index(variabletag2_list[i])])
			
	densityplot_MERGER_SAT_NORM(xlistall, ylistall, xnames, ynames, redshiftnames, xlims, ylims, "densityplot_z_merger_sat_" + plottag+"NORM.pdf", 40e0,plottag)

def flyby(traj):
	for i in range(len(traj))[4:]:
		if traj[i-3] < traj[i-2] < traj[i-1]:
			if traj[i-3] < traj[i-4] and traj[i-1] < traj[i]:
#				if traj[i-2] < d_flyby:
				return i
	return -1

def findflyby(filelist0A, filelist0B):
	plotlist_r = []
	for l in range(len(z_redshift_list)):
		dist_satA = plotread(filelist0A, l, 'dist_sat')[2]
		dist_satB = plotread(filelist0B, l, 'dist_sat')[2]
		dist_sat = np.concatenate((dist_satA, dist_satB))
		
		plotlist_r.append(dist_sat)
	
	plotlist_r_T = np.transpose(np.array(plotlist_r))
	
	f = open('flybyAB_d_flyby.dat','w')
	f2 = open('NEGflybyAB_d_flyby.dat','w')
	flyby_index_list = []
	flyby_neg_index_list = []
	for j in range(len(plotlist_r_T)):
		i = flyby(plotlist_r_T[j])
		if i > 0:
			f.write(str(j) + '\t' + str(i) + '\n')
			flyby_index_list.append(j)
		else:
			f2.write(str(j) +'\n')
			flyby_neg_index_list.append(j)
	f.close()
	f2.close()

def individualsat(filelist0):
	#plot individual satellite radius and cos angle vs redshift (try to see flybys)

	fig = plt.figure()
	fig.set_size_inches(10, 10, forward=True)
	gs = gridspec.GridSpec(2,2)
	gs.update(left=0.1, right=0.9, hspace=0.2, wspace=0.25)

	gridnum = 2
	
	axlist = []
	for j in range(2):
		for i in range(2):
			axlist.append(plt.subplot(gs[j, i]))

	cmap = mpl.cm.rainbow
	plotlist_r = []
	plotlist_cos = []
	for l in range(len(z_redshift_list)):
		sat_cosangle_pos = plotread(filelist0, l, 'sat_cosangle_pos')[2] #rp
		dist_sat = plotread(filelist0, l, 'dist_sat')[2]		
		plotlist_r.append(dist_sat)
		plotlist_cos.append(sat_cosangle_pos)

	plotlist_r = np.transpose(np.array(plotlist_r))
	plotlist_cos = np.transpose(np.array(plotlist_cos))

	flyby_index_list = loadtxt('flybyAB_d_flyby.dat', dtype = int, usecols = (0,))
	flyby_neg_index_list = loadtxt('NEGflybyAB_d_flyby.dat', dtype = int)

	num_traj = int(50/gridnum) #number of "trajectories", radius vs. time

	for j in range(gridnum):
		for i in range(num_traj):
			axlist[j].plot(z_redshift_list, plotlist_r[flyby_index_list[:2500:50]][j * num_traj + i], color=cmap(i / float(20)))
			axlist[j].set_ylim(0, 1)
			axlist[j].set_xlabel("z",fontsize=18)
			axlist[j].set_ylabel('$d_{\mathrm{sat}}$',fontsize=18)
			
	for j in range(gridnum):
		for i in range(num_traj):
			axlist[j+2].plot(z_redshift_list, plotlist_cos[flyby_index_list[:2500:50]][j * num_traj + i], color=cmap(i / float(20)))
			axlist[j+2].set_xlabel("z",fontsize=18)
			axlist[j+2].set_ylabel('$\cos\\theta$',fontsize=18)
	fig.savefig('./paper/individual_radius_cos.pdf', dpi=500)
	plt.close()

def minradiusdensity(filelist0A, filelist0B):
	#plot individual satellite radius and cos angle vs redshift (try to see flybys)
	fig = plt.figure(figsize=(6,6))
	ax = fig.add_subplot(111)

	flyby_position_list = loadtxt('flybyAB_d_flyby.dat', dtype = int, usecols = (1, ))
	radius_bar = []
	for l in range(len(z_redshift_list)):
		zind = np.where(np.array(flyby_position_list) == l)
		
		radius_bar.append(len(zind[0]))
	
	labels = ['flyby event', 'cos position angle']
	colors = ['red', 'blue']
	ax.bar(np.array(z_redshift_list), radius_bar, color=colors[0], width = 0.1, label = labels[0])

	fig.savefig('density_minradius.pdf', dpi=500)
	plt.close()
	
def merger_sat_quiver_plot(filelist0, plotname, zoom):
#	z_plotlist = [0.0, 0.0175, 0.5179, 0.9489, 1.5025, 2.5511]
#	z_plotlist = [0.0, 0.0175, 0.0696, 0.1245, 0.3065,0.5179]
#	z_plotlist = [0.0, 0.0696,0.2427, 0.5957, 1.1542,1.5025, 1.7663, 2.0572,2.5511]
	z_plotlist = [0.0, 0.5179, 0.9489, 1.5025, 1.9078,2.5511]
	redshiftnames = []
	l_list = []
	for l in range(27): # over z
		if z_redshift_list[l] in z_plotlist: # over z
			redshiftnames.append('z = ' + str(round(z_redshift_list[l],2)))
			l_list.append(l)
			
	fig = plt.figure()
	fig.set_size_inches(14,10, forward=True)
	gs = gridspec.GridSpec(2,3)
	gs.update(left=0.1, right=0.9, hspace=0.2, wspace=0.2)
	
	axlist = []
	for j in range(2):
		for i in range(3):
			axlist.append(plt.subplot(gs[j, i]))
	axlist[1].set(aspect = 1, title = 'projected satellite density')
	
	for j in range(len(redshiftnames)):
		sat_cosangle_pos = plotread(filelist0, l_list[j], 'sat_cosangle_pos')[2] #rp
		sat_sinangle_pos = plotread(filelist0, l_list[j], 'sat_sinangle_pos')[2] #rp
		pv_relhost_cosangle = plotread(filelist0, l_list[j], 'pv_relhost_cosangle')[2] #vr
		v_relhost_cosangle = plotread(filelist0, l_list[j], 'sat_cosangle_vel_relhost')[2] #vp
		
		dist_sat = plotread(filelist0, l_list[j], 'dist_sat')[2]		
		vsat = plotread(filelist0, l_list[j], 'vsat')[2]
		logmass_sat = plotread(filelist0, l_list[j], 'logmass_sat')[2]
		
		theta = np.random.rand(len(sat_cosangle_pos))*2e0*pi

		x_sat = dist_sat * sat_cosangle_pos
		y_sat = dist_sat * sat_sinangle_pos * cos(theta)
		
#		xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
#		ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
		
#		xmax = max(x_sat)
#		xmin = min(x_sat)

		if zoom == 0:
			xmax = 3.5e0
			xmin = -3.5e0
			ymax = 3.5e0
			ymin = -3.5e0
			nbins = 26
			vmin0 = 1
			vmax0 = 4
			
		if zoom == 1:
			xmax = 2e0
			xmin = -2e0
			ymax = 2e0
			ymin = -2e0
			nbins = 26
			vmax0 = 3
			vmin0 = 1
			
		deltax = (xmax-xmin)/(2e0 * nbins)
		deltay = (ymax-ymin)/(2e0 * nbins)
		
		x0 = np.arange(xmin, xmax + deltax, deltax)
		y0 = np.arange(ymax, ymin - deltay, -deltay)
			
		x = np.arange(xmin + deltax/2e0, xmax, deltax) # plot grids center
		y = np.arange(ymax - deltay/2e0, ymin, -deltay)
		X, Y = np.meshgrid(x,y)
		
		density_mtr = []
		for i in range(nbins*2): #x direction grid
			density_mtr.append([])
			for k in range(nbins*2): #y direction grid
				index = np.where( np.logical_and(np.logical_and(x_sat > x0[k], x_sat < x0[k + 1]), np.logical_and(y_sat < y0[i], y_sat > y0[i + 1])))
				density_mtr[i].append(np.log10(len(index[0])))

		axlist[j].set_xlabel(redshiftnames[j])
		
		im = axlist[j].imshow(np.array(density_mtr), extent=[xmin, xmax, ymin, ymax], cmap='jet', aspect = 'equal', interpolation='nearest', vmin=vmin0, vmax=vmax0)
		axlist[j].contour(X, Y, np.array(density_mtr), colors='black')
		axlist[j].set(aspect = 1, title = 'z = ' + str(round(z_redshift_list[l_list[j]],2)))
	
		axlist[j].grid(True, linestyle= '-', which= 'major',color= '0.75')
		axlist[j].grid(True, which='both')
		axlist[j].set_xlabel('x')
		axlist[j].set_ylabel('y')
							
	fig.colorbar(im, format='$%.2f$',cax = fig.add_axes([0.91, 0.2, 0.02, 0.6]))
	if zoom == 1:	
		fig.savefig(os.path.join(os.getcwd(), "paper", plotname+'density_zoomin.pdf'), dpi=500)
	elif zoom == 0:
		fig.savefig(os.path.join(os.getcwd(), "paper", plotname+'density_zoomout.pdf'), dpi=500)
	plt.close()

def merger_sat_quiver_plot_ID(filelist0, plotname, ID_tag, zoom):
#	z_plotlist = [0.0, 0.0696, 0.1822, 0.2427, 0.5957, 1.1542,1.5025, 1.7663, 2.0572]
#	z_plotlist = [0.0, 0.5179, 0.9489, 1.5025, 1.9078, 2.5511]
#	z_plotlist = [0.0, 0.0175, 0.5179, 0.9489, 1.5025, 2.5511]
	z_plotlist = [0.0, 0.0175, 0.0696,0.1245,0.1822, 0.2427,0.3734,0.5179]

	l_list = []
	for l in range(27): # over z
		if z_redshift_list[l] in z_plotlist: # over z
			l_list.append(l)
			
	fig = plt.figure()
	fig.set_size_inches(18, 12, forward=True)
	gs = gridspec.GridSpec(2,4)
	gs.update(left=0.1, right=0.9, hspace=0.15, wspace=0.2)
	
	axlist = []
	for j in range(2):
		for i in range(4):
			axlist.append(plt.subplot(gs[j, i]))

	if zoom == 1:
		xmax = 1.
		xmin = -1.
		ymax = 1.
		ymin = -1.
					
		nbins = 26
		vmax0 = 2000
		
	if zoom == 0:
		xmax = 2.
		xmin = -2.
		ymax = 2.
		ymin = -2
		nbins = 26
		vmax0 = 2000
				
	plotlist_r = []
	plotlist_cos = []
	for j in range(len(l_list)):
		axlist[j].set(aspect = 1, title = 'z = ' + str(round(z_redshift_list[l_list[j]],2)))
		sat_cosangle_pos = plotread(filelist0, l_list[j], 'sat_cosangle_pos')[2]
		sat_sinangle_pos = plotread(filelist0, l_list[j], 'sat_sinangle_pos')[2] #rp
		dist_sat = plotread(filelist0, l_list[j], 'dist_sat')[2]
		if j == 0:
			if ID_tag == 0: # ID_tag = 0, face forward, pie slice, = 1, face backward from partner, pie slice, = 2, extrema, =3, extrema
				index = np.where(sat_cosangle_pos > 0.8)
			###################
			elif ID_tag == 1:
				index = np.where(sat_cosangle_pos < -0.8)
			###################	
			elif ID_tag == 2:
				index = np.where(np.logical_and(sat_cosangle_pos > 0.8, dist_sat > 0.4 ))
			###################
			elif ID_tag == 3:
				index = np.where(np.logical_and(sat_cosangle_pos < -0.8, dist_sat > 0.4 )) 
				###################
			index = index[0]
			
			theta = np.random.rand(len(index))*2e0*pi

		x_sat_tot = dist_sat[index] * sat_cosangle_pos[index]
		y_sat_tot = dist_sat[index] * sat_sinangle_pos[index] * cos(theta)
				
		deltax = (xmax-xmin)/(2e0 * nbins)
		deltay = (ymax-ymin)/(2e0 * nbins)
		
		x0 = np.arange(xmin, xmax + deltax, deltax)
		y0 = np.arange(ymax, ymin - deltay, -deltay)
			
		x = np.arange(xmin + deltax/2e0, xmax, deltax) # plot grids center
		y = np.arange(ymax - deltay/2e0, ymin, -deltay)
		X, Y = np.meshgrid(x,y)
				
		density_mtr = []
		for i in range(nbins * 2): #x direction grid
			density_mtr.append([])
			for k in range(nbins * 2): #y direction grid
				index0 = np.where( np.logical_and(np.logical_and(x_sat_tot > x0[k], x_sat_tot < x0[k + 1]), np.logical_and(y_sat_tot < y0[i], y_sat_tot > y0[i + 1])))
				density_mtr[i].append(len(index0[0])/len(index))
		axlist[j].contour(X, Y, np.array(density_mtr), colors='black')
		im = axlist[j].imshow(np.array(density_mtr), extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Blues, aspect = 'equal', interpolation='nearest')		
		axlist[j].set_xlabel('x')
		axlist[j].set_ylabel('y')

		if zoom == 1:
			if ID_tag == 0: # ID_tag = 0, face forward, pie slice, = 1, face backward from partner, pie slice, = 2, extrema, =3, extrema
				fig.savefig(plotname+'_face_partner_zoomin_smallz.pdf', dpi=500)	
			elif ID_tag == 1:
				fig.savefig(plotname+'_opp_partner_zoomin_smallz.pdf', dpi=500)
			elif ID_tag == 2:
				fig.savefig(plotname+'_extr_face_partner_zoomin_smallz.pdf', dpi=500)
			elif ID_tag == 3:
				fig.savefig(plotname+'_extr_opp_partner_zoomin_smallz.pdf', dpi=500)
		elif zoom == 0:
			if ID_tag == 0: # ID_tag = 0, face forward, pie slice, = 1, face backward from partner, pie slice, = 2, extrema, =3, extrema
				fig.savefig(plotname+'_face_partner_zoomout_smallz.pdf', dpi=500)	
			elif ID_tag == 1:
				fig.savefig(plotname+'_opp_partner_zoomout_smallz.pdf', dpi=500)
			elif ID_tag == 2:
				fig.savefig(plotname+'_extr_face_partner_zoomout_smallz.pdf', dpi=500)
			elif ID_tag == 3:
				fig.savefig(plotname+'_extr_opp_partner_zoomout_smallz.pdf', dpi=500)		
	plt.close()

def merger_sat_quiver_plot_ID_flyby(filelist0A, filelist0B, plotname, zoom, ID_tag):
#	z_plotlist = [0.0, 0.5179, 0.9489, 1.5025, 1.9078, 2.5511]
	z_plotlist = [0.0, 0.0175, 0.1245,0.1822,0.3734,0.5179]
	l_list = []
	for l in range(27): # over z
		if z_redshift_list[l] in z_plotlist: # over z
			l_list.append(l)
			
	fig = plt.figure()
	fig.set_size_inches(14, 12, forward=True)
	gs = gridspec.GridSpec(2,3)
	gs.update(left=0.1, right=0.9, hspace=0.15, wspace=0.2)
	
	axlist = []
	for j in range(2):
		for i in range(3):
			axlist.append(plt.subplot(gs[j, i]))

	xmax = 1.
	xmin = -1.
	ymax = 1.
	ymin = -1.
				
	nbins = 26
	vmax0 = 200
		
	plotlist_r = []
	plotlist_cos = []
	plotlist_sin = []
	
	for l in range(len(z_redshift_list)):
		sat_cosangle_pos_A = plotread(filelist0A, l, 'sat_cosangle_pos')[2]
		sat_cosangle_pos_B = plotread(filelist0B, l, 'sat_cosangle_pos')[2]
		sat_cosangle_pos = np.concatenate((sat_cosangle_pos_A, sat_cosangle_pos_B))	
		if l == 0:
			if ID_tag == 0: # ID_tag = 0, face forward, pie slice, = 1, face backward from partner, pie slice, = 2, extrema, =3, extrema
				index = np.where(sat_cosangle_pos > 0.8)
			###################
			elif ID_tag == 1:
				index = np.where(sat_cosangle_pos < -0.8)
			###################	
			elif ID_tag == 2:
				index = np.where(np.logical_and(sat_cosangle_pos > 0.8, dist_sat > 0.4 ))
			###################
			elif ID_tag == 3:
				index = np.where(np.logical_and(sat_cosangle_pos < -0.8, dist_sat > 0.4 )) 
				###################
			index = index[0]
			
		sat_sinangle_pos_A = plotread(filelist0A, l, 'sat_sinangle_pos')[2]
		sat_sinangle_pos_B = plotread(filelist0B, l, 'sat_sinangle_pos')[2]
		sat_sinangle_pos = np.concatenate((sat_sinangle_pos_A, sat_sinangle_pos_B))
		
		dist_satA = plotread(filelist0A, l, 'dist_sat')[2]
		dist_satB = plotread(filelist0B, l, 'dist_sat')[2]
		dist_sat = np.concatenate((dist_satA, dist_satB))
		
		plotlist_r.append(dist_sat)
		plotlist_cos.append(sat_cosangle_pos)
		plotlist_sin.append(sat_sinangle_pos)
		
	plotlist_r_T = np.transpose(np.array(plotlist_r))
	
	flyby_index_list0 = loadtxt('flybyAB_d_flyby.dat', dtype = int, usecols = (0,))
#	flyby_neg_index_list = loadtxt('NEGflybyAB_d_flyby.dat', dtype = int)
	
	indflyby = np.where( np.in1d(flyby_index_list0, index) )
	print(len(flyby_index_list0), len(index), len(indflyby[0]))
	flyby_index_list = flyby_index_list0[indflyby[0]]
	
	theta = np.random.rand( len(flyby_index_list) ) * 2e0 * pi
	
	for j in range(len(l_list)):
		axlist[j].set(aspect = 1, title = 'z = ' + str(round(z_redshift_list[l_list[j]],2)))
		dist_sat = plotlist_r[l_list[j]]
		sat_cosangle_pos = plotlist_cos[l_list[j]]
		sat_sinangle_pos = plotlist_sin[l_list[j]]
		x_sat_tot = dist_sat[flyby_index_list] * sat_cosangle_pos[flyby_index_list]
		y_sat_tot = dist_sat[flyby_index_list] * sat_sinangle_pos[flyby_index_list] * cos(theta)					
		deltax = (xmax-xmin)/(2e0 * nbins)
		deltay = (ymax-ymin)/(2e0 * nbins)
		
		x0 = np.arange(xmin, xmax + deltax, deltax)
		y0 = np.arange(ymax, ymin - deltay, -deltay)
			
		x = np.arange(xmin + deltax/2e0, xmax, deltax) # plot grids center
		y = np.arange(ymax - deltay/2e0, ymin, -deltay)
		X, Y = np.meshgrid(x,y)
				
		density_mtr = []
		for i in range(nbins * 2): #x direction grid
			density_mtr.append([])
			for k in range(nbins * 2): #y direction grid
				index0 = np.where( np.logical_and(np.logical_and(x_sat_tot > x0[k], x_sat_tot < x0[k + 1]), np.logical_and(y_sat_tot < y0[i], y_sat_tot > y0[i + 1])))
				density_mtr[i].append(np.log10(len(index0[0])))
				
		im = axlist[j].imshow(np.array(density_mtr), extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Blues, aspect = 'equal', interpolation='nearest')		

		axlist[j].set_xlabel('x')
		axlist[j].set_ylabel('y')
		axlist[j].contour(X, Y, np.array(density_mtr), colors='black', alpha = 0.5)
		axlist[j].grid(True, linestyle= '-', which= 'major',color= '0.75')
		axlist[j].grid(True, which='both')
		fig.colorbar(im, format='$%.2f$',cax = fig.add_axes([0.91, 0.2, 0.02, 0.6]))	
		if ID_tag == 0: # ID_tag = 0, face forward, pie slice, = 1, face backward from partner, pie slice, = 2, extrema, =3, extrema
			fig.savefig(plotname+'_face_partner_zoomin_smallz.pdf', dpi=500)	
		elif ID_tag == 1:
			fig.savefig(plotname+'_opp_partner_zoomin_smallz.pdf', dpi=500)
		elif ID_tag == 2:
			fig.savefig(plotname+'_extr_face_partner_zoomin_smallz.pdf', dpi=500)
		elif ID_tag == 3:
			fig.savefig(plotname+'_extr_opp_partner_zoomin_smallz.pdf', dpi=500)
	plt.close()

def signal_flyby(filelist0A, filelist0B, plotname, slicetag):
	
	plotlist_r = []
	plotlist_cos = []
					
	for l in range(len(z_redshift_list)):
		sat_cosangle_pos_A = plotread(filelist0A, l, 'sat_cosangle_pos')[2]
		sat_cosangle_pos_B = plotread(filelist0B, l, 'sat_cosangle_pos')[2]
		sat_cosangle_pos = np.concatenate((sat_cosangle_pos_A, sat_cosangle_pos_B))

		dist_satA = plotread(filelist0A, l, 'dist_sat')[2]
		dist_satB = plotread(filelist0B, l, 'dist_sat')[2]
		dist_sat = np.concatenate((dist_satA, dist_satB))
		
		plotlist_r.append(dist_sat)
		plotlist_cos.append(sat_cosangle_pos)
		
	plotlist_r_T = np.transpose(np.array(plotlist_r))
	
	flyby_index_list = loadtxt('flybyAB_d_flyby.dat', dtype = int, usecols = (0,))
	flyby_neg_index_list = loadtxt('NEGflybyAB_d_flyby.dat', dtype = int)

	if slicetag == 0:
		dist_sat = plotlist_r[0][flyby_index_list]
		sat_cosangle_pos = plotlist_cos[0][flyby_index_list]
		dist_sat_2 = plotlist_r[0][flyby_neg_index_list]
		sat_cosangle_pos_2 = plotlist_cos[0][flyby_neg_index_list]
		dist_sat_3 = plotlist_r[0] # signal
		sat_cosangle_pos_3 = plotlist_cos[0]

	else:
		dist_sat = plotlist_r[19][flyby_index_list]
		sat_cosangle_pos = plotlist_cos[19][flyby_index_list]
		dist_sat_2 = plotlist_r[19][flyby_neg_index_list]
		sat_cosangle_pos_2 = plotlist_cos[19][flyby_neg_index_list]
		dist_sat_3 = plotlist_r[19] # signal
		sat_cosangle_pos_3 = plotlist_cos[19]
		
	Xlist = sat_cosangle_pos
	Ylist = dist_sat
	Xlist2 = sat_cosangle_pos_2
	Ylist2 = dist_sat_2
	Xlist3 = sat_cosangle_pos_3
	Ylist3 = dist_sat_3
	
	print(len(Xlist), len(Xlist2), len(Xlist3))
#	sys.exit()
	plottitles = ['(a) Flybys' , '(b) Non-flybys', '(c) All Satellites']
	if slicetag == 0:
		densityplot_3panels(Xlist, Ylist, Xlist2, Ylist2, Xlist3, Ylist3, '$\cos\\theta$', '$d_{\mathrm{sat}}/d_{\mathrm{sep0}}$', plotname, plottitles, [-1e0,1e0], [0,.5e0], [0.4, 0.8], 2.8, 40)
	else:
		densityplot_3panels(Xlist, Ylist, Xlist2, Ylist2, Xlist3, Ylist3, '$\cos\\theta$', '$d_{\mathrm{sat}}/d_{\mathrm{sep0}}$', plotname, plottitles, [-1e0,1e0], [0,1.2e0], [0.33, 1.5], 1.3, 40)


def signal_inoutstreaming(filelist0A, filelist0B, plotname, slicetag):
	
	sat_cosangle_pos_A = plotread(filelist0A, 0, 'sat_cosangle_pos')[2]
	sat_cosangle_pos_B = plotread(filelist0B, 0, 'sat_cosangle_pos')[2]
	sat_cosangle_pos0 = np.concatenate((sat_cosangle_pos_A, sat_cosangle_pos_B))

	dist_satA = plotread(filelist0A, 0, 'dist_sat')[2]
	dist_satB = plotread(filelist0B, 0, 'dist_sat')[2]
	dist_sat0 = np.concatenate((dist_satA, dist_satB))
	
	pv_relhost_cosangle_A = plotread(filelist0A, 0, 'pv_relhost_cosangle')[2] #vr
	pv_relhost_cosangle_B = plotread(filelist0B, 0, 'pv_relhost_cosangle')[2] #vr
	pv_relhost_cosangle = np.concatenate((pv_relhost_cosangle_A, pv_relhost_cosangle_B))	

	OUT_index_list = np.where(pv_relhost_cosangle > 0)
	IN_index_list = np.where(pv_relhost_cosangle < 0)
	dist_sat = dist_sat0[OUT_index_list[0]]
	sat_cosangle_pos = sat_cosangle_pos0[OUT_index_list[0]]
	dist_sat_2 = dist_sat0[IN_index_list[0]]
	sat_cosangle_pos_2 = sat_cosangle_pos0[IN_index_list[0]]
	dist_sat_3 = dist_sat0 # signal
	sat_cosangle_pos_3 = sat_cosangle_pos0

	Xlist = sat_cosangle_pos
	Ylist = dist_sat
	Xlist2 = sat_cosangle_pos_2
	Ylist2 = dist_sat_2
	Xlist3 = sat_cosangle_pos_3
	Ylist3 = dist_sat_3
	
	plottitles = ['(a) Streaming out ' , '(b) Streaming in', '(c) All Satellites']
	
	if slicetag == 0:
		densityplot_3panels(Xlist, Ylist, Xlist2, Ylist2, Xlist3, Ylist3, '$\cos\\theta$', '$d_{\mathrm{sat}}/d_{\mathrm{sep0}}$', plotname, plottitles, [-1e0,1e0], [0,.5e0], [0.4, 0.74], 2.8, 40)
	else:
		densityplot_3panels(Xlist, Ylist, Xlist2, Ylist2, Xlist3, Ylist3, '$\cos\\theta$', '$d_{\mathrm{sat}}/d_{\mathrm{sep0}}$', plotname, plottitles, [-1e0,1e0], [0,1.2e0], [0.33, 1.5], 1.3, 40)

		
def merger_sat_density_plot_flyby(filelist0A, filelist0B, plottag, facingtag, flybytag):
	z_plotlist = [0,0.0175,0.1245,0.3734,1.0492,2.5511]

	for l in range(len(z_redshift_list)):
		if l == 0:
			sat_cosangle_posA = plotread(filelist0A, l, 'sat_cosangle_pos')[2]
			sat_cosangle_posB = plotread(filelist0B, l, 'sat_cosangle_pos')[2]	
			sat_cosangle_pos = np.concatenate((sat_cosangle_posA, sat_cosangle_posB))
			if facingtag == 'facing':
				index0 = np.where(sat_cosangle_pos > 0.8)
			else:
				index0 = np.where(sat_cosangle_pos < - 0.8)
						
			index0 = index0[0]	
	
	if flybytag == 'flyby':
		flyby_index_list0 = loadtxt('flybyAB_d_flyby.dat', dtype = int, usecols = (0,))
		indflyby = np.where( np.in1d(flyby_index_list0, index0) )
		flyby_index_list = flyby_index_list0[indflyby[0]]
	else:
		flyby_neg_index_list0 = loadtxt('NEGflybyAB_d_flyby.dat', dtype = int)
		indNEGflyby = np.where( np.in1d(flyby_neg_index_list0, index0) )
		flyby_neg_index_list = flyby_neg_index_list0[indNEGflyby[0]]
	
#	flyby_neg_index_list = flyby_neg_index_list0
	
	xlistall = []
	ylistall = []
	xlims = []
	ylims = []
	xnames = []
	ynames = []
	counter = 0
	redshiftnames = []
	
	for l in range(27): # over z
		if z_redshift_list[l] in z_plotlist: # over z
			xlistall.append([])
			ylistall.append([])
			counter += 1
			redshiftnames.append('z = ' + str(round(z_redshift_list[l],2)))

			for i in range(len(variabletag1_list)): # over variables to plot
				pmax1, pmin1, var_list1_A = plotread(filelist0A, l, variabletag1_list[i]) # massiveornot = 0 for massive, 1 for less massive
				pmax2, pmin2, var_list2_A = plotread(filelist0A, l, variabletag2_list[i])
				
				pmax1, pmin1, var_list1_B = plotread(filelist0B, l, variabletag1_list[i]) # massiveornot = 0 for massive, 1 for less massive
				pmax2, pmin2, var_list2_B = plotread(filelist0B, l, variabletag2_list[i])
				
				var_list1 = np.concatenate((var_list1_A, var_list1_B))
				var_list2 = np.concatenate((var_list2_A, var_list2_B))
				if flybytag == 'flyby':
					xlistall[counter - 1].append(var_list1[flyby_index_list])
					ylistall[counter - 1].append(var_list2[flyby_index_list])
				else:
					xlistall[counter - 1].append(var_list1[flyby_neg_index_list])
					ylistall[counter - 1].append(var_list2[flyby_neg_index_list])				
				if l == 0:
					xlims.append([pmin1, pmax1])
					ylims.append([pmin2, pmax2])
					xnames.append(ylabel[varname_list.index(variabletag1_list[i])])
					ynames.append(ylabel[varname_list.index(variabletag2_list[i])])
						
	densityplot_MERGER_SAT_NORM(xlistall, ylistall, xnames, ynames, redshiftnames, xlims, ylims, "densityplot_z_merger_sat_" + plottag+"NORM_"+ flybytag + '_' + facingtag+"2.pdf", 40e0,plottag)

def merger_sat_density_plot_flyby2(filelist0A, filelist0B, plottag, facingtag, flybytag):
	z_plotlist = [0,0.0175,0.1245,0.3734,1.0492,2.5511]

	for l in range(len(z_redshift_list)):
		if l == 0:
			sat_cosangle_posA = plotread(filelist0A, l, 'sat_cosangle_pos')[2]
			sat_cosangle_posB = plotread(filelist0B, l, 'sat_cosangle_pos')[2]	
			sat_cosangle_pos = np.concatenate((sat_cosangle_posA, sat_cosangle_posB))
			if facingtag == 'facing':
				index0 = np.where(sat_cosangle_pos > 0.8)
			else:
				index0 = np.where(sat_cosangle_pos < - 0.8)
						
			index0 = index0[0]
	
	flyby_index_list0 = loadtxt('flybyAB_d_flyby.dat', dtype = int, usecols = (0,))
	indflyby = np.where( np.in1d(flyby_index_list0, index0) )
	flyby_index_list = flyby_index_list0[indflyby[0]]
	flyby_neg_index_list0 = loadtxt('NEGflybyAB_d_flyby.dat', dtype = int)
	indNEGflyby = np.where( np.in1d(flyby_neg_index_list0, index0) )
	flyby_neg_index_list = flyby_neg_index_list0[indNEGflyby[0]]
	
	xlistall = []
	ylistall = []
	xlims = []
	ylims = []
	xnames = []
	ynames = []
	counter = 0
	redshiftnames = []
	
	for l in range(27): # over z
		if z_redshift_list[l] in z_plotlist: # over z
			xlistall.append([])
			ylistall.append([])
			counter += 1
			redshiftnames.append('z = ' + str(round(z_redshift_list[l],2)))

			for i in range(2): # over variables to plot
				pmax1, pmin1, var_list1_A = plotread(filelist0A, l, variabletag1_list[0]) # massiveornot = 0 for massive, 1 for less massive
				pmax2, pmin2, var_list2_A = plotread(filelist0A, l, variabletag2_list[0])
				
				pmax1, pmin1, var_list1_B = plotread(filelist0B, l, variabletag1_list[0]) # massiveornot = 0 for massive, 1 for less massive
				pmax2, pmin2, var_list2_B = plotread(filelist0B, l, variabletag2_list[0])
				
				var_list1 = np.concatenate((var_list1_A, var_list1_B))
				var_list2 = np.concatenate((var_list2_A, var_list2_B))

				if i == 0:
					ylistall[counter - 1].append(var_list1[flyby_index_list])
					xlistall[counter - 1].append(var_list2[flyby_index_list])
				else:
					ylistall[counter - 1].append(var_list1[flyby_neg_index_list])
					xlistall[counter - 1].append(var_list2[flyby_neg_index_list])				
				if l == 0:
					ylims.append([pmin1, pmax1])
					xlims.append([pmin2, pmax2])
					ynames.append(ylabel[varname_list.index(variabletag1_list[i])])
					xnames.append(ylabel[varname_list.index(variabletag2_list[i])])
						
	densityplot_MERGER_SAT_NORM2(xlistall, ylistall, xnames, ynames, redshiftnames, xlims, ylims, "densityplot_z_merger_sat_" + plottag+"NORM_"+ flybytag + '_' + facingtag+"3.pdf", 40e0,plottag)
	
def merger_sat_density_plot_position(filelist0, plottag):
	z_plotlist = [0.0, 0.5179, 1.5025,2.5511]
#	z_plotlist = [0.0, 0.0175, 0.0696,0.1245]
#	z_plotlist = [0.1245,0.3734,0.6776,1.0492]
#	l_list = []
#	for l in range(len(z_redshift_list)): # over z
#		if z_redshift_list[l] in z_plotlist: # over z
#			l_list.append(l)
			
	plotlist_r = []

	xlistall = []
	ylistall = []
	xlims = []
	ylims = []
	xnames = []
	ynames = []
	counter = 0
	redshiftnames = []
	for j in range(len(z_redshift_list)):
		if j == 0:
			dist_sat = plotread(filelist0, j, 'dist_sat')[2]
	
			sat_cosangle_pos = plotread(filelist0, j, 'sat_cosangle_pos')[2]		
			index0 = np.where(sat_cosangle_pos > 0.8)
			index0 = index0[0]
			
		if z_redshift_list[j] in z_plotlist:
			l = j
			xlistall.append([])
			ylistall.append([])
			
			redshiftnames.append('z = ' + str(round(z_redshift_list[l],1)))
			counter += 1
			for i in range(len(variabletag1_list)): # over variables to plot
				pmax1, pmin1, var_list1 = plotread(filelist0, l, variabletag1_list[i]) # massiveornot = 0 for massive, 1 for less massive
				pmax2, pmin2, var_list2 = plotread(filelist0, l, variabletag2_list[i])
			
				xlistall[counter - 1].append(var_list1[index0])
				ylistall[counter - 1].append(var_list2[index0])
			
				if counter == 1:
					xlims.append([pmin1, pmax1])
					ylims.append([pmin2, pmax2])
					xnames.append(ylabel[varname_list.index(variabletag1_list[i])])
					ynames.append(ylabel[varname_list.index(variabletag2_list[i])])
			
	print(xnames, xlistall[0])
#	densityplot_MERGER_SAT_NORM(xlistall, ylistall, xnames, ynames, redshiftnames, xlims, ylims, "densityplot_z_merger_sat_" + plottag+"NORM_facing_opp_midz.pdf", 40e0,plottag)
	densityplot_MERGER_SAT_NORM(xlistall, ylistall, xnames, ynames, redshiftnames, xlims, ylims, os.path.join(os.getcwd(), "paper", "densityplot_z_merger_sat_" + plottag+"NORM_facing.pdf"), 40e0,plottag)			

def vel_facing_vs_opp(filelist0A, filelist0B):
	fig = plt.figure()
	ax = fig.add_subplot(111)

	flyby_index_list0 = loadtxt('flybyAB_d_flyby.dat', dtype = int, usecols = (0,))
	NEGflyby_index_list0 = loadtxt('NEGflybyAB_d_flyby.dat', dtype = int)

#	sat_cosangle_pos_A = plotread(filelist0A, 0, 'sat_cosangle_pos')[2]
#	sat_cosangle_pos_B = plotread(filelist0B, 0, 'sat_cosangle_pos')[2]	
	sat_cosangle_pos_A = plotread(filelist0A, len(z_redshift_list)-1, 'sat_cosangle_pos')[2]
	sat_cosangle_pos_B = plotread(filelist0B, len(z_redshift_list)-1, 'sat_cosangle_pos')[2]
	sat_cosangle_pos = np.concatenate((sat_cosangle_pos_A, sat_cosangle_pos_B))
	index0_facing = np.where(sat_cosangle_pos > 0)
	index_facing = index0_facing[0]
	index0_facing_opp = np.where(sat_cosangle_pos < 0)
	index_facing_opp = index0_facing_opp[0]

	vsat_avg_facing_list = []
	vsat_avg_facing_opp_list = []
	vsat_avg_facing_list2 = []
	vsat_avg_facing_opp_list2 = []	
	
	for l in range(len(z_redshift_list)):
		vsat_A = plotread(filelist0A, l, 'vsat')[2]
		vsat_B = plotread(filelist0B, l, 'vsat')[2]
		vsat = np.concatenate((vsat_A, vsat_B))

		indflyby_fac = np.where( np.in1d(flyby_index_list0, index_facing) )
		flyby_index_list_fac = flyby_index_list0[indflyby_fac[0]]
		
		indflyby_opp = np.where( np.in1d(flyby_index_list0, index_facing_opp) )
		flyby_index_list_opp = flyby_index_list0[indflyby_opp[0]]

		indNEGflyby_fac = np.where( np.in1d(NEGflyby_index_list0, index_facing) )
		NEGflyby_index_list_fac = NEGflyby_index_list0[indNEGflyby_fac[0]]
		
		indNEGflyby_opp = np.where( np.in1d(NEGflyby_index_list0, index_facing_opp) )
		NEGflyby_index_list_opp = NEGflyby_index_list0[indNEGflyby_opp[0]]
						
		vsat_avg_facing = np.mean(vsat[flyby_index_list_fac])
		vsat_avg_facing_opp = np.mean(vsat[flyby_index_list_opp])															
		vsat_avg_facing_list.append(vsat_avg_facing)
		vsat_avg_facing_opp_list.append(vsat_avg_facing_opp)

		vsat_avg_facing2 = np.mean(vsat[NEGflyby_index_list_fac])
		vsat_avg_facing_opp2 = np.mean(vsat[NEGflyby_index_list_opp])															
		vsat_avg_facing_list2.append(vsat_avg_facing2)
		vsat_avg_facing_opp_list2.append(vsat_avg_facing_opp2)
		
	ax.plot(z_redshift_list, vsat_avg_facing_list, color='red', label = 'facing partner, flyby')
	ax.plot(z_redshift_list, vsat_avg_facing_opp_list, color='blue', label = 'facing away from partner, flyby')
	ax.plot(z_redshift_list, vsat_avg_facing_list2, color='red',  linestyle = '--', label = 'facing partner, nonflyby')
	ax.plot(z_redshift_list, vsat_avg_facing_opp_list2, color='blue',  linestyle = '--', label = 'facing away from partner, nonflyby')	
	ax.set_xlabel("z",fontsize=20)
	ax.set_xlim(0, z_redshift_list[-1])
	ax.set_ylabel("$<v_{\mathrm{sat}}>$", fontsize=20)	
	ax.grid(True, linestyle= '-', which= 'major',color= '0.75')
	ax.grid(True, linestyle= '-', which= 'minor',color= '0.75')
	ax.grid(True, which='both')
	ax.minorticks_on()
	legend = ax.legend(loc="upper right", shadow=True, fontsize = 15)
	fig.savefig('vel_facing_vs_opp.pdf', dpi=500)
	plt.close()	

def vel_facing_vs_opp2(filelist0, hosttag):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	flyby_index_list0 = loadtxt('flybyAB_d_flyby.dat', dtype = int, usecols = (0,))
	NEGflyby_index_list0 = loadtxt('NEGflybyAB_d_flyby.dat', dtype = int)

	sat_cosangle_pos = plotread(filelist0, len(z_redshift_list)-1, 'sat_cosangle_pos')[2]
	
	index0_facing = np.where(sat_cosangle_pos > 0)
	index_facing = index0_facing[0]
	index0_facing_opp = np.where(sat_cosangle_pos < 0)
	index_facing_opp = index0_facing_opp[0]

	vsat_avg_facing_list = []
	vsat_avg_facing_opp_list = []
	vsat_avg_facing_list2 = []
	vsat_avg_facing_opp_list2 = []	

	cos_avg_facing_list = []
	cos_avg_facing_opp_list = []
	cos_avg_facing_list2 = []
	cos_avg_facing_opp_list2 = []		
	for l in range(len(z_redshift_list)):
		vsat = plotread(filelist0, l, 'vsat')[2]
		sat_cosangle_pos = plotread(filelist0, l, 'sat_cosangle_vel_relhost')[2]
		
		indflyby_fac = np.where( np.in1d(flyby_index_list0, index_facing) )
		flyby_index_list_fac = flyby_index_list0[indflyby_fac[0]]
		
		indflyby_opp = np.where( np.in1d(flyby_index_list0, index_facing_opp) )
		flyby_index_list_opp = flyby_index_list0[indflyby_opp[0]]

		indNEGflyby_fac = np.where( np.in1d(NEGflyby_index_list0, index_facing) )
		NEGflyby_index_list_fac = NEGflyby_index_list0[indNEGflyby_fac[0]]
		
		indNEGflyby_opp = np.where( np.in1d(NEGflyby_index_list0, index_facing_opp) )
		NEGflyby_index_list_opp = NEGflyby_index_list0[indNEGflyby_opp[0]]
						
#		vsat_avg_facing = np.mean(vsat[flyby_index_list_fac])
#		vsat_avg_facing_opp = np.mean(vsat[flyby_index_list_opp])															
#		vsat_avg_facing_list.append(vsat_avg_facing)
#		vsat_avg_facing_opp_list.append(vsat_avg_facing_opp)

#		vsat_avg_facing2 = np.mean(vsat[NEGflyby_index_list_fac])
#		vsat_avg_facing_opp2 = np.mean(vsat[NEGflyby_index_list_opp])															
#		vsat_avg_facing_list2.append(vsat_avg_facing2)
#		vsat_avg_facing_opp_list2.append(vsat_avg_facing_opp2)

		cos_avg_facing = np.mean(sat_cosangle_pos[flyby_index_list_fac])
		cos_avg_facing_opp = np.mean(sat_cosangle_pos[flyby_index_list_opp])															
		cos_avg_facing_list.append(cos_avg_facing)
		cos_avg_facing_opp_list.append(cos_avg_facing_opp)
		
		cos_avg_facing2 = np.mean(sat_cosangle_pos[NEGflyby_index_list_fac])
		cos_avg_facing_opp2 = np.mean(sat_cosangle_pos[NEGflyby_index_list_opp])															
		cos_avg_facing_list2.append(cos_avg_facing2)
		cos_avg_facing_opp_list2.append(cos_avg_facing_opp2)
#				
#	ax.plot(z_redshift_list, vsat_avg_facing_list, color='red', label = 'facing partner, flyby')
#	ax.plot(z_redshift_list, vsat_avg_facing_opp_list, color='blue', label = 'facing away from partner, flyby')
#	ax.plot(z_redshift_list, vsat_avg_facing_list2, color='red',  linestyle = '--', label = 'facing partner, nonflyby')
#	ax.plot(z_redshift_list, vsat_avg_facing_opp_list2, color='blue',  linestyle = '--', label = 'facing away from partner, nonflyby')
	
	ax.plot(z_redshift_list, cos_avg_facing_list, color='orange', label = 'facing partner, flyby')
	ax.plot(z_redshift_list, cos_avg_facing_opp_list, color='green', label = 'facing away from partner, flyby')
	ax.plot(z_redshift_list, cos_avg_facing_list2, color='orange',  linestyle = '--', label = 'facing partner, nonflyby')
	ax.plot(z_redshift_list, cos_avg_facing_opp_list2, color='green',  linestyle = '--', label = 'facing away from partner, nonflyby')	
		
	ax.set_xlabel("z",fontsize=20)
	ax.set_xlim(0, z_redshift_list[-1])
	ax.set_ylabel("$<v_{\mathrm{sat}}>$", fontsize=20)	
	ax.grid(True, linestyle= '-', which= 'major',color= '0.75')
	ax.grid(True, linestyle= '-', which= 'minor',color= '0.75')
	ax.grid(True, which='both')
	ax.minorticks_on()
	legend = ax.legend(loc="upper right", shadow=True, fontsize = 15)
	fig.savefig('cosvel_facing_vs_opp' + hosttag + '.pdf', dpi=500)
	plt.close()	

def vel_facing_vs_opp3(filelist0A, filelist0B):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	color_list = ['red', 'blue', 'orange', 'green']
	for j, hosttag in enumerate(['A','B']):
		filelist0 = [filelist0A, filelist0B][j]
		flyby_index_list0 = loadtxt('flybyAB_d_flyby.dat', dtype = int, usecols = (0,))
		NEGflyby_index_list0 = loadtxt('NEGflybyAB_d_flyby.dat', dtype = int)

		sat_cosangle_pos = plotread(filelist0, len(z_redshift_list)-1, 'sat_cosangle_pos')[2]
	
		index0_facing = np.where(sat_cosangle_pos > 0)
		index_facing = index0_facing[0]
		index0_facing_opp = np.where(sat_cosangle_pos < 0)
		index_facing_opp = index0_facing_opp[0]

		vsat_avg_facing_list = []
		vsat_avg_facing_opp_list = []
		vsat_avg_facing_list2 = []
		vsat_avg_facing_opp_list2 = []	

		cos_avg_facing_list = []
		cos_avg_facing_opp_list = []
		cos_avg_facing_list2 = []
		cos_avg_facing_opp_list2 = []		
		for l in range(len(z_redshift_list)):
			vsat = plotread(filelist0, l, 'vsat')[2]
			sat_cosangle_pos = plotread(filelist0, l, 'sat_cosangle_vel_relhost')[2]
		
			indflyby_fac = np.where( np.in1d(flyby_index_list0, index_facing) )
			flyby_index_list_fac = flyby_index_list0[indflyby_fac[0]]
		
			indflyby_opp = np.where( np.in1d(flyby_index_list0, index_facing_opp) )
			flyby_index_list_opp = flyby_index_list0[indflyby_opp[0]]

			indNEGflyby_fac = np.where( np.in1d(NEGflyby_index_list0, index_facing) )
			NEGflyby_index_list_fac = NEGflyby_index_list0[indNEGflyby_fac[0]]
		
			indNEGflyby_opp = np.where( np.in1d(NEGflyby_index_list0, index_facing_opp) )
			NEGflyby_index_list_opp = NEGflyby_index_list0[indNEGflyby_opp[0]]
						
			vsat_avg_facing = np.mean(vsat[flyby_index_list_fac])
			vsat_avg_facing_opp = np.mean(vsat[flyby_index_list_opp])															
			vsat_avg_facing_list.append(vsat_avg_facing)
			vsat_avg_facing_opp_list.append(vsat_avg_facing_opp)

			vsat_avg_facing2 = np.mean(vsat[NEGflyby_index_list_fac])
			vsat_avg_facing_opp2 = np.mean(vsat[NEGflyby_index_list_opp])															
			vsat_avg_facing_list2.append(vsat_avg_facing2)
			vsat_avg_facing_opp_list2.append(vsat_avg_facing_opp2)

#			cos_avg_facing = np.mean(sat_cosangle_pos[flyby_index_list_fac])
#			cos_avg_facing_opp = np.mean(sat_cosangle_pos[flyby_index_list_opp])															
#			cos_avg_facing_list.append(cos_avg_facing)
#			cos_avg_facing_opp_list.append(cos_avg_facing_opp)
#		
#			cos_avg_facing2 = np.mean(sat_cosangle_pos[NEGflyby_index_list_fac])
#			cos_avg_facing_opp2 = np.mean(sat_cosangle_pos[NEGflyby_index_list_opp])															
#			cos_avg_facing_list2.append(cos_avg_facing2)
#			cos_avg_facing_opp_list2.append(cos_avg_facing_opp2)
#		
		ax.plot(z_redshift_list, vsat_avg_facing_list, color=color_list[j * 2], label = 'facing partner at $z=2.6$, flyby, ' + hosttag)
		ax.plot(z_redshift_list, vsat_avg_facing_opp_list, color=color_list[j*2+1], label = 'facing away from partner at $z=2.6$, flyby, ' + hosttag)
		ax.plot(z_redshift_list, vsat_avg_facing_list2, color=color_list[j * 2],  linestyle = '--', label = 'facing partner at $z=2.6$, nonflyby, ' + hosttag)
		ax.plot(z_redshift_list, vsat_avg_facing_opp_list2, color=color_list[j*2+1],  linestyle = '--', label = 'facing away from partner at $z=2.6$, nonflyby, ' + hosttag)
#		ax.plot(z_redshift_list, cos_avg_facing_list, color=color_list[j * 2], label = 'facing partner at $z=2.6$, flyby, ' + hosttag)
#		ax.plot(z_redshift_list, cos_avg_facing_opp_list, color=color_list[j*2+1], label = 'facing away from partner at $z=2.6$, flyby, ' + hosttag)
#		ax.plot(z_redshift_list, cos_avg_facing_list2, color=color_list[j * 2],  linestyle = '--', label = 'facing partner at $z=2.6$, nonflyby, ' + hosttag)
#		ax.plot(z_redshift_list, cos_avg_facing_opp_list2,color=color_list[j * 2+1],linestyle = '--', label = 'facing away from partner at $z=2.6$, nonflyby, ' + hosttag)
	ax.set_xlabel("z",fontsize=20)
	ax.set_xlim(0, z_redshift_list[-1])
#	ax.set_ylabel('$\cos\\theta_{\mathrm{vel}}$', fontsize=20)
	ax.set_ylabel("$<v_{\mathrm{sat}}>$", fontsize=20)	
	ax.grid(True, linestyle= '-', which= 'major',color= '0.75')
	ax.grid(True, linestyle= '-', which= 'minor',color= '0.75')
	ax.grid(True, which='both')
	ax.minorticks_on()
#	legend = ax.legend(loc="center right", shadow=True, fontsize = 10)
#	fig.savefig('cosvel_facing_vs_opp_AB.pdf', dpi=500)
	
	legend = ax.legend(loc="upper right", shadow=True, fontsize = 10)
	fig.savefig('vel_facing_vs_opp_AB.pdf', dpi=500)
	plt.close()	

def diff_flyby_nonflyby(filelist0):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	color_list = ['red', 'blue']
	
	flyby_index_list = loadtxt('flybyA.dat', dtype = int, usecols = (0,))
	flyby_neg_index_list = loadtxt('NEGflybyA.dat', dtype = int)	
	
	plot_mass_flyby = []
	plot_mass_nonflyby = []
	print(z_redshift_list)
	for l in range(len(z_redshift_list)):
		vsat = plotread(filelist0, l, 'vsat')[2]
		logmass_sat = plotread(filelist0, l, 'logmass_sat')[2]
#		print(filelist0[l])
#		print(logmass_sat[flyby_index_list][:10])
		plot_mass_flyby.append( mean(logmass_sat[flyby_index_list]))
		plot_mass_nonflyby.append( mean(logmass_sat[flyby_neg_index_list]))
#		print('mass', mean(logmass_sat[flyby_index_list]), mean(logmass_sat[flyby_neg_index_list]))
#		print('v', mean(vsat[flyby_index_list]), mean(vsat[flyby_neg_index_list]))
	ax.plot(z_redshift_list, plot_mass_flyby, color=color_list[0], label = 'flyby')
	ax.plot(z_redshift_list, plot_mass_nonflyby, color=color_list[1], label = 'nonflyby')
	ax.set_xlabel("z",fontsize=20)
	ax.set_xlim(0, z_redshift_list[-1])
	ax.set_ylabel("log mass", fontsize=20)	
	ax.grid(True, linestyle= '-', which= 'major',color= '0.75')
	ax.grid(True, linestyle= '-', which= 'minor',color= '0.75')
	ax.grid(True, which='both')
	ax.minorticks_on()	
	
	legend = ax.legend(loc="upper right", shadow=True, fontsize = 10)
	fig.savefig('log_mass_flyby.pdf', dpi=500)
	plt.close()	

def overlappedABdensity(filelist0A, filelist0B):
	z_plotlist = [0.0, 0.5179, 0.9489, 1.5025, 1.9078,2.5511]
	redshiftnames = []
	l_list = []
	for l in range(27): # over z
		if z_redshift_list[l] in z_plotlist: # over z
			redshiftnames.append('z = ' + str(round(z_redshift_list[l],2)))
			l_list.append(l)
			
	fig = plt.figure()
	fig.set_size_inches(19,10, forward=True)
	gs = gridspec.GridSpec(2,3)
	gs.update(left=0.1, right=0.9, hspace=0.15, wspace=0.15)
	
	axlist = []
	for j in range(2):
		for i in range(3):
			axlist.append(plt.subplot(gs[j, i]))
	axlist[1].set(aspect = 1, title = 'projected satellite density')
	
	A_avg_logmass_list = []
	B_avg_logmass_list = []
	for j in range(len(redshiftnames)):
		#### calculate average host mass per redshift slice
#		logm_max, logm_min, log_mass = plotread_host(filelist0_host, l_list[j], 'logmass')
#		A_logmass = log_mass[::2]
#		B_logmass = log_mass[1::2]
#		major_logmass = []
#		minor_logmass = []
#		for i in range(len(A_logmass)):
#			if A_logmass[i] > B_logmass[i]:
#				major_logmass.append(A_logmass[i])
#				minor_logmass.append(B_logmass[i])
#			else:
#				major_logmass.append(B_logmass[i])
#				minor_logmass.append(A_logmass[i])
				
#		A_avg_logmass_1z = mean(major_logmass)
#		B_avg_logmass_1z = mean(minor_logmass)
#		A_avg_logmass_list.append(A_avg_logmass_1z)
#		B_avg_logmass_list.append(B_avg_logmass_1z)
		
		sat_cosangle_pos_A = plotread(filelist0A, l_list[j], 'sat_cosangle_pos')[2] #rp
		sat_sinangle_pos_A = plotread(filelist0A, l_list[j], 'sat_sinangle_pos')[2] #rp
		sat_cosangle_pos_B = plotread(filelist0B, l_list[j], 'sat_cosangle_pos')[2] #rp
		sat_sinangle_pos_B = plotread(filelist0B, l_list[j], 'sat_sinangle_pos')[2] #rp		
		pv_relhost_cosangle_A = plotread(filelist0A, l_list[j], 'pv_relhost_cosangle')[2] #vr
		v_relhost_cosangle_A = plotread(filelist0A, l_list[j], 'sat_cosangle_vel_relhost')[2] #vp
		pv_relhost_cosangle_B = plotread(filelist0B, l_list[j], 'pv_relhost_cosangle')[2] #vr
		v_relhost_cosangle_B = plotread(filelist0B, l_list[j], 'sat_cosangle_vel_relhost')[2]
		
		dist_sat_A = plotread(filelist0A, l_list[j], 'dist_sat')[2]		
		dist_sat_B = plotread(filelist0B, l_list[j], 'dist_sat')[2]		
				
		theta_A = np.random.rand(len(sat_cosangle_pos_A))*2e0*pi

		theta_B = np.random.rand(len(sat_cosangle_pos_B))*2e0*pi
		
		x_satA = dist_sat_A * sat_cosangle_pos_A
		y_satA = dist_sat_A * sat_sinangle_pos_A * cos(theta_A)

		x_satB = np.ones(len(sat_cosangle_pos_B)) - dist_sat_B * sat_cosangle_pos_B
		y_satB = dist_sat_B * sat_sinangle_pos_B * cos(theta_B)
		
		x_sat = np.concatenate((x_satA, x_satB))
		y_sat = np.concatenate((y_satA, y_satB))

		xmax = 4e0
		xmin = -4e0
		ymax = 4e0
		ymin = -4e0
						
#		xmax = 2e0
#		xmin = -1e0
#		ymax = 1e0
#		ymin = -1e0
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
		
#		density_mtrA = []
#		for i in range(nbins*2): #x direction grid
#			density_mtrA.append([])
#			for k in range(nbins*2): #y direction grid
#				index = np.where( np.logical_and(np.logical_and(x_satA > x0[k], x_satA < x0[k + 1]), np.logical_and(y_satA < y0[i], y_satA > y0[i + 1])))
#				density_mtrA[i].append(np.log10(len(index[0])))

#		im = axlist[j].imshow(np.array(density_mtrA), extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Reds, aspect = 'equal', interpolation='nearest', vmin=vmin0, vmax=vmax0, alpha = 0.7)
#		axlist[j].contour(X, Y, np.array(density_mtrA), colors='red', alpha = 0.3,levels =np.arange(.6,4.,.4))
#		
#		density_mtrB = []
#		for i in range(nbins*2): #x direction grid
#			density_mtrB.append([])
#			for k in range(nbins*2): #y direction grid
#				index = np.where( np.logical_and(np.logical_and(x_satB > x0[k], x_satB < x0[k + 1]), np.logical_and(y_satB < y0[i], y_satB > y0[i + 1])))
#				density_mtrB[i].append(np.log10(len(index[0])))

#		im = axlist[j].imshow(np.array(density_mtrB), extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Blues, aspect = 'equal', interpolation='nearest', vmin=vmin0, vmax=vmax0, alpha = 0.7)
#		axlist[j].contour(X, Y, np.array(density_mtrB), colors='blue', alpha = 0.3,levels =np.arange(.6,4.,.4))


		density_mtr = []
		for i in range(nbins*2): #x direction grid
			density_mtr.append([])
			for k in range(nbins*2): #y direction grid
				index = np.where( np.logical_and(np.logical_and(x_sat > x0[k], x_sat < x0[k + 1]), np.logical_and(y_sat < y0[i], y_sat > y0[i + 1])))
				density_mtr[i].append(np.log10(len(index[0])))

#		im = axlist[j].imshow(np.array(density_mtr), extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Blues, aspect = 'equal', interpolation='nearest', vmin=vmin0, vmax=vmax0)
		im = axlist[j].imshow(np.array(density_mtr), extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Blues, aspect = 'equal', interpolation='nearest')		
		axlist[j].contour(X, Y, np.array(density_mtr), colors='black')

		axlist[j].set(aspect = 1, title = 'z = ' + str(round(z_redshift_list[l_list[j]],2)))
		plt.colorbar(im, format='$%.2f$', ax = axlist[j])
		axlist[j].grid(True, linestyle= '-', which= 'major',color= '0.75')
		axlist[j].grid(True, which='both')
		axlist[j].set_xlabel('x')
		axlist[j].set_ylabel('y')
		
#	t = np.arange(vmin0,vmax0,.4)
#	fig.colorbar(im, format='$%.2f$',cax = fig.add_axes([0.91, 0.2, 0.02, 0.6]))
#	fig.savefig('ABdensity_overlapped_2color_zoomout.pdf', dpi=500)
	fig.savefig('ABdensity_overlapped_zoomin.pdf', dpi=500)
	fig.savefig('ABdensity_overlapped_zoomout.pdf', dpi=500)			
	plt.close()
		
#	DAT = np.column_stack((np.array(z_plotlist), A_avg_logmass_list, B_avg_logmass_list))
#	np.savetxt('host_avg_logmass_vs_z.dat', DAT, delimiter=" ")


def velvector(potential, i, j, deltax, deltay):
	V_x = potential[i][j+1]
	V2_x = potential[i][j-1]
	u_x = (V2_x - V_x)/(2e0*deltax)

	V_y = potential[i-1][j]
	V2_y = potential[i+1][j]
	u_y = (V2_y - V_y)/(2e0*deltay)
	return u_x, u_y
	
def AB_potential(normalize = 1):
#	deltax = 5e-2
#	deltay = 5e-2

#	xmin = -4e0
#	xmax = 4e0
#	ymin = -4e0
#	ymax = 4e0
#	countourres = 3e-1

	deltax = 4e-2
	deltay = 4e-2	
	xmin = -2e0
	xmax = 2e0
	ymin = -2e0
	ymax = 2e0
	xylim = [xmin, xmax, ymin, ymax]
	countourres = 2e-2
	
	x_mid = np.arange(xmin + deltax/2e0, xmax, deltax)
	y_mid = np.arange(ymax - deltay/2e0, ymin, - deltay)
	Xmid, Ymid = np.meshgrid(x_mid, y_mid)	 
	Xmid2, Ymid2 = np.meshgrid(x_mid-1e0, y_mid)
	
	z, log10mA, log10mB = loadtxt('host_avg_logmass_vs_z.dat', unpack = 1)
	mA = (10e0**(log10mA[0]))/(10e0**(log10mB[0]))
	mB = 1e0
	print(mA,mB) # 2.66517033464 = major host, 1.0 = minor host
	
	A_pot = - mA / sqrt(Xmid ** 2e0 + Ymid ** 2e0)
	B_pot = - mB / sqrt(Xmid2 ** 2e0 + Ymid2 ** 2e0)
	
	fig = plt.figure()
	ax = fig.add_subplot(111)

	ax.set_xlabel("x",fontsize=20)
	ax.set_ylabel("y", fontsize=20)	
	
	AB_pot = A_pot + B_pot
	im = ax.imshow(AB_pot, extent=[xmin, xmax, ymin, ymax], cmap=plt.cm.Blues, aspect = 'equal', interpolation='nearest', vmin=-9., vmax=-.4)
	U = []
	V = []
	
	print(len(AB_pot),len(AB_pot[0]))
	for i in range(1,len(AB_pot)-1):
		U.append([])
		V.append([])
		for j in range(1,len(AB_pot[0])-1):
#			print(i,j)
			[ux, uy] = velvector(AB_pot, i, j, deltax, deltay) 
			speed = sqrt(ux ** 2e0 + uy ** 2e0)
			if normalize == 0:
				U[i-1].append(ux)
				V[i-1].append(uy)
			if normalize == 1:
				U[i-1].append(ux/log2(speed))
				V[i-1].append(uy/log2(speed))		
		U[i-1] = np.array(U[i-1])
		V[i-1] = np.array(V[i-1])
	U = np.array(U)
	V = np.array(V)
	plt.quiver(Xmid[::5,::5], Ymid[::5,::5], U[::5,::5], V[::5,::5], pivot='mid', scale= 50)
	
	mpl.rcParams['contour.negative_linestyle'] = 'solid'
	ax.contour(Xmid, Ymid, AB_pot, colors='black', vmin=-1, vmax=-.42, alpha = 0.3, levels =-(1e0/(np.arange(1e-2, 2e0, countourres))))
	
	if normalize == 0:
		fig.savefig('ABpotential_quiver.jpg', dpi=500)
	else:
		fig.savefig('ABpotential_quiver_speed.jpg', dpi=500)

	plt.close()		

	return AB_pot, U, V, Xmid, Ymid, deltax, deltay, xylim
	
def theoretical_signal():
	AB_pot, U, V, Xmid, Ymid, deltax, deltay, xylim = AB_potential(normalize = 0)

	[xmin, xmax, ymin, ymax] = xylim

	N = 100000
	r3d = np.random.rand(N)
	r = 2e0 * np.cbrt(r3d) #### !!!! this is to ensure the volumn is uniformly sampled, not just radius
	phi = np.random.rand(N) * 2e0 * pi
	costheta = np.random.rand(N) * 2e0 - 1e0

	xarr = r * costheta
	yarr = r * sqrt(1e0 - costheta ** 2e0) * cos(phi)
		
	Nint = 1000
	h = 1e-3

	x_edges = np.arange(xmin, xmax + deltax, deltax)
	y_edges = np.arange(ymax, ymin - deltay, - deltay)
#	print(np.shape(x_edges))
#	print(np.shape(y_edges))	
#	print(np.shape(U))	
	for i in range(Nint):
		xarr2 = []
		yarr2 = []	
		for j in range(len(xarr)):
#			if xarr[j] ** 2e0 + yarr[j] ** 2e0 < 0.1:
#				continue
#			elif (xarr[j]-1e0) ** 2e0 + yarr[j] ** 2e0 < 0.1:
#				continue				
#			else:
			for k in range(1, len(x_edges) - 2):
#				print(k, j)
				if xarr[j] < x_edges[k+1] and xarr[j] > x_edges[k]:
					for l in range(1, len(y_edges) - 2):
						if yarr[j] > y_edges[l+1] and yarr[j] < y_edges[l]:
							newx = xarr[j] + h * U[l-1, k-1]
							newy = yarr[j] + h * V[l-1, k-1]
							if np.sign(newx) != np.sign(xarr[j]) or np.sign(newy) != np.sign(yarr[j]) or np.sign(abs(newx-1e0)) != np.sign(abs(xarr[j]-1e0)):
								continue
#								print('append1')
#								xarr2.append(0)
#								yarr2.append(0)
							else:
#								print('append2')
								xarr2.append(newx)
								yarr2.append(newy)
		xarr = xarr2
		yarr = yarr2

#		indexX = np.where(np.logical_or(np.abs(xarr) < 1e-1, np.abs(xarr - 1e0) < 1e-1) )
#		indexX = indexX[0]
#		indexY = np.where(np.abs(yarr[indexX]) < 1e-1)
#		indexY = indexY[0]
		
#		yarr = yarr[indexX[indexY]]
#		xarr = xarr[indexX[indexY]]
		
#		print(yarr2[:100])
				
#		for j in range(1, len(x_edges) - 2):
#			index0 = np.where(np.logical_and(xarr[indexX[indexY]] > x_edges[j], xarr[indexX[indexY]] < x_edges[j+1]))
#			index0 = index0[0]
#			for k in range(1, len(y_edges) - 2):
#				index1 = np.where(np.logical_and(xarr[indexX[indexY]][index0] < y_edges[k], yarr2[index0] > y_edges[k+1]))
##				print(xarr[index0[index1[0]]], h * U[j-1, k-1])
#				xarr[indexX[indexY]] = h * U[j-1, k-1]
#				yarr[] = h * V[j-1, k-1]	
##				print(xarr[index0[index1[0]]])

		if int(i % 10) == 0:	
			H, xedges, yedges = np.histogram2d(yarr, xarr, bins=(y_edges[::-1], x_edges))
			fig = plt.figure(figsize=(8,8))
			ax = fig.add_subplot(111)

			im = ax.imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', aspect = 'equal', vmin=0, vmax=50)
			plt.colorbar(im, format='$%.2f$',cax = fig.add_axes([0.91, 0.2, 0.02, 0.6]))
			plt.savefig('theoretical_signal_time=' + str(i) + '.jpg', transparent=True, dpi=500,bbox_inches='tight')
			plt.close()	
			
if __name__ == "__main__":
###################################################################################################################################################################
	colorlist = []
	for name, hex in matplotlib.colors.cnames.items():
		colorlist.append(name)

	filelist0_host = []
	for file_1z in os.listdir('./'):
		if file_1z.endswith("_mmax.dat"):
			filelist0_host.append(file_1z)
#		if file_1z.endswith("_mmax_prognum.dat"):
#			filelist0_prog.append(file_1z)
	filelist0_host.sort(key = natural_keys)
			
	filelist0A = []
	filelist0B = []
	for file_1z in os.listdir('./'):
		if file_1z.endswith("samesizesample_hostA.dat"):
#		if file_1z.endswith("merger_info_hostA.dat"):
			filelist0A.append(file_1z)
		elif file_1z.endswith("samesizesample_hostB.dat"):
#		elif file_1z.endswith("merger_info_hostB.dat"):
			filelist0B.append(file_1z)
	filelist0A.sort(key = natural_keys)
	filelist0B.sort(key = natural_keys)
	
#	for j, filelist0 in enumerate([filelist0A, filelist0B]):
#		for l in range(len(filelist0)):
##			print(filelist0[l])
#			pmax, pmin, m_list = plotread(filelist0, l, 'logmass_sat')
##			print(len(m_list))
##			print(( 10e0 ** m_list[-100000] )/particlemass)
#			toosmallhalo = (( 10e0 ** m_list )/particlemass < 20)
#			idofbigenoughhalos = np.where(toosmallhalo == 0)
##			print(len(idofbigenoughhalos[0]))
#			if l != 0:
#				toosmallhalo = toosmallhalo_old + toosmallhalo
#			toosmallhalo_old = toosmallhalo

#			print(10e0**min(m_list[idofbigenoughhalos[0]])/particlemass)
#			
#		idofbigenoughhalos = np.where(toosmallhalo == 0)
#		
#		print(len(idofbigenoughhalos[0]))
#		if j == 0:
#			np.savetxt('idofbigenoughhalosA.dat' , idofbigenoughhalos[0], delimiter=" ")		
#		else:
#			np.savetxt('idofbigenoughhalosB.dat' , idofbigenoughhalos[0], delimiter=" ")		
#	sys.exit()
######################################################################################	
	z_redshift_list = []
	for l in range(len(filelist0A)):
		z_redshift = float(filelist0A[l][3:9])
		z_redshift_list.append(z_redshift)

	z_redshift_list_linear = np.arange(0, z_redshift_list[-1], .3)
	
	ylabel = ['$\log(M_{sat}$)', '$\cos\\theta$', '$\cos\\theta_{\mathrm{vel}}$', '$\cos\\theta_{vel relcent}$', '$\cos\\theta_{\mathrm{rv}}$', '$\cos\\theta_{pv rel cent}$', '$d_{\mathrm{sat}}/d_{\mathrm{sep0}}$', '$v_{\mathrm{sat}}$']
	varname_list = ['logmass_sat', 'sat_cosangle_pos', 'sat_cosangle_vel_relhost', 'sat_cosangle_vel_relcent', 'pv_relhost_cosangle', 'pv_relcent_cosangle', 'dist_sat', 'vsat']
	
	ylabel2 = ['$\cos\\theta$', '$\cos\\theta_{\mathrm{rv}}$']
	varname_list2 = ['sat_cosangle_pos', 'pv_relhost_cosangle']

#	ylabel2 = [ '$\cos\\theta_{\mathrm{rv}}$']
#	varname_list2 = [ 'pv_relhost_cosangle']
		
################################################################################	
#MERGER_LINES

#	for l in range(len(filelist0B)):
#		print('redshift:', z_redshift_list[l])
#		
#		pmax, pmin, m_list = plotread(filelist0B,l,'logmass_sat')
#		print(amin(m_list), len(m_list))
#		
#		sys.exit()
#		if l == 0:
#			idn = np.where(m_list == amin(m_list))
#			ID_list = plotread(filelist0A,l,'ID')
#			idid = ID_list[idn[0]][0]
#			ID = int(ID_list[idn[0]][0])
#			print('ID=', ID)
#		
#		if l > 0:
#			descID_list = plotread(filelist0A,l,'descID')
#			
#			print('descID=', int(descID_list[idn[0]][0]))
#			
#			iddid = np.where(descID_list == ID)
#			
#			ID_list = plotread(filelist0A,l,'ID')
##			print('other haloes sharing the same descID:', int(ID_list[iddid[0]]))
#			ID = int(ID_list[idn[0]][0])
#			print('ID=', ID)		
#		
#		print('mass = ', 10** m_list[idn[0]][0], "particle N: ",  10 ** m_list[idn[0]][0]/2.6e6)
#		
#		print('\n')
		
################################################################################
################################################################################		
#	for vartag in varname_list:
#	for vartag in ['logmass_sat']:
#		merger_sat_analysis_plot(filelist0A, vartag, 'hostA')
#		merger_sat_analysis_plot(filelist0B, vartag, 'hostB')
#	sys.exit()
##	
#	for vartag in varname_list2: #slice plots
###		merger_sat_analysis_plot3d(vartag, 'A')
###		merger_sat_analysis_plot3d(vartag, 'B')
#		merger_sat_analysis_plot3d(vartag, 'AB')	

#	sys.exit()
#	variabletag1_list = ['sat_cosangle_pos', 'pv_relhost_cosangle']
#	variabletag2_list = ['dist_sat', 'dist_sat']

	variabletag1_list = ['pv_relhost_cosangle']
	variabletag2_list = ['dist_sat']
	
#	merger_sat_density_plot(filelist0A, filelist0B, 'AB')
#	sys.exit()

#	merger_sat_density_plot(filelist0A, 'hostA')
#	merger_sat_density_plot(filelist0B, 'hostB')

#	merger_sat_quiver_plot(filelist0A, 'quiverplotA', 1)
#	merger_sat_quiver_plot(filelist0B, 'quiverplotB', 1)
#	merger_sat_quiver_plot(filelist0A, 'quiverplotA', 0)
#	merger_sat_quiver_plot(filelist0B, 'quiverplotB', 0)
#	sys.exit()
#	merger_sat_quiver_plot_ID(filelist0A, 'quiverplotA', 0, 1)
#	merger_sat_quiver_plot_ID(filelist0B, 'quiverplotB', 0, 1)
#	merger_sat_quiver_plot_ID(filelist0A, 'quiverplotA', 1, 1)
#	merger_sat_quiver_plot_ID(filelist0B, 'quiverplotB', 1, 1)

#	merger_sat_quiver_plot_ID(filelist0A, 'quiverplotA', 2, 1)
#	merger_sat_quiver_plot_ID(filelist0B, 'quiverplotB', 2, 1)
#	merger_sat_quiver_plot_ID(filelist0A, 'quiverplotA', 3, 1)
#	merger_sat_quiver_plot_ID(filelist0B, 'quiverplotB', 3, 1)

##	merger_sat_quiver_plot_ID(filelist0A, 'quiverplotA', 0, 0)
##	merger_sat_quiver_plot_ID(filelist0B, 'quiverplotB', 0, 0)
##	merger_sat_quiver_plot_ID(filelist0A, 'quiverplotA', 1, 0)
##	merger_sat_quiver_plot_ID(filelist0B, 'quiverplotB', 1, 0)
##	merger_sat_quiver_plot_ID(filelist0A, 'quiverplotA', 2, 0)
##	merger_sat_quiver_plot_ID(filelist0B, 'quiverplotB', 2, 0)
##	merger_sat_quiver_plot_ID(filelist0A, 'quiverplotA', 3, 0)
##	merger_sat_quiver_plot_ID(filelist0B, 'quiverplotB', 3, 0)

########## flyby

	variabletag1_list = ['dist_sat', 'dist_sat', 'dist_sat']
	variabletag2_list = ['sat_cosangle_pos', 'sat_cosangle_vel_relhost','vsat']

#	merger_sat_quiver_plot_ID_flyby(filelist0A, filelist0B, 'flyby', 0, 0)
#	merger_sat_quiver_plot_ID_flyby(filelist0A, filelist0B, 'flyby', 1, 1)
#	sys.exit()

#	merger_sat_density_plot_flyby(filelist0A, filelist0B, 'hostA+B', 'facing', 'flyby')
#	merger_sat_density_plot_flyby(filelist0A, filelist0B, 'hostA+B', 'facing_opp', 'flyby')
#	merger_sat_density_plot_flyby(filelist0A, filelist0B, 'hostA+B', 'facing', 'nonflyby')
#	merger_sat_density_plot_flyby(filelist0A, filelist0B, 'hostA+B', 'facing_opp', 'nonflyby')

	findflyby(filelist0A, filelist0B)
	
#	signal_flyby(filelist0A, filelist0B, os.path.join(os.getcwd(), "paper", "rtheta_densityplot_signal_flyby_3panels_z=1.5_2.pdf"),19)
	signal_flyby(filelist0A, filelist0B, os.path.join(os.getcwd(), "paper", "rtheta_densityplot_signal_flyby_3panels_z=0_2.pdf"),0)		
#	merger_sat_density_plot_flyby2(filelist0A, filelist0B, 'hostA+B', 'facing', 'flyby+nonflyby')
#	sys.exit()
	


#	signal_inoutstreaming(filelist0A, filelist0B, os.path.join(os.getcwd(), "paper", "rtheta_densityplot_signal_streaming_3panels_z=0_2.pdf"),0)	
	
#	individualsat(filelist0A)
#	minradiusdensity(filelist0A, filelist0B)

#	merger_sat_density_plot_position(filelist0A, 'hostA')
#	merger_sat_density_plot_position(filelist0B, 'hostB')

#	findflyby2(filelist0A, 'A')
#	findflyby2(filelist0B, 'B')	

#	vel_facing_vs_opp(filelist0A, filelist0B)
#	vel_facing_vs_opp2(filelist0A, 'A')
#	vel_facing_vs_opp2(filelist0B, 'B')	
#	vel_facing_vs_opp3(filelist0A, filelist0B)		

#	diff_flyby_nonflyby(filelist0A)
#	overlappedABdensity(filelist0A, filelist0B)
#	AB_potential(normalize = 0)
#	theoretical_signal()
