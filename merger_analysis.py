'''
This code studies the merger tree history of the LGB of today z = 0
'''
import matplotlib as mpl
from numpy import loadtxt, histogram
import numpy as np
import sys
from constant import *
import os
from shutil import copyfile
from helpers import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
from numpy import amax, amin, median, std
from matplotlib.ticker import NullFormatter, MaxNLocator, AutoMinorLocator

def RockstarLGBpath(Rockstarfilelist, filenum):
	filename_high_res = Rockstarfilelist[filenum]
	print(filename_high_res)
	a = float(filename_high_res[6:12]) #scale
	z_redshift = (1.- a)/a #redshift
	
	# for each redshift there is a folder to store data in
	datfoldr = "/snapshot_for_z=%.4f" % z_redshift
	os.chdir(os.getcwd() + datfoldr)
	LGB_path = os.getcwd() + "/LGB_z=%.4f.dat" % z_redshift
	
	return LGB_path, z_redshift
	
def merger_analysis():
	LGBfilename0, z_redshift = RockstarLGBpath(Rockstarfilelist, 0)
	IDold, descIDold, m, x, y, z, vx, vy, vz, rvir, b, c = col_reader6(LGBfilename0, IDcol, descIDcol, masscol, xcol, ycol, zcol, vxcol, vycol, vzcol, Rvircol, bcol, ccol)
	DAT = np.column_stack((IDold, m, vx, vy, vz, x, y, z, rvir, b, c))
	np.savetxt('z_=%.4f_mmax.dat' % z_redshift, DAT, delimiter=" ")
	os.chdir("..")
	deadend_num_list = []
	
	for l in range(1, 30): #for each redshift outside z = 0
#		mostmassive_prog_num = []
		deadend = 0
		LGBfilename, z_redshift = RockstarLGBpath(Rockstarfilelist, l)
		IDnew, descIDnew, m, x, y, z, vx, vy, vz, rvir, b, c = col_reader6(LGBfilename, IDcol, descIDcol, masscol, xcol, ycol, zcol, vxcol, vycol, vzcol, Rvircol, bcol, ccol)
		
		mmax_list_1z = [] #max mass list for 1 redshift value
		mmax_index_list_1z = []
		
		for i in range(len(IDold)):
			index = np.where( np.in1d( descIDnew, IDold[i] ) )
#			
			if len(index[0]) !=	0:
				mmax = np.amax(m[index[0]]) #pick out biggest predecessor mass
				mmax_list_1z.append(mmax) #append to biggest predecessor mass list
				for j in range(len(index[0])):
					if int(mmax) == int(m[index[0][j]]):
						mmax_index = index[0][j]
				mmax_index_list_1z.append(mmax_index)
#				mostmassive_prog_num.append(len(index[0]))
			else:
				deadend += 1
		deadend_num_list.append(deadend)
#		np.savetxt('z=%.4f_mmax_prognum.dat' % z_redshift, mostmassive_prog_num, delimiter=" ")
		
#		print(descIDnew[mmax_index], IDold[i])
		IDnew_list_1z = IDnew[mmax_index_list_1z]
		vxmax_list_1z = vx[mmax_index_list_1z]
		vymax_list_1z = vy[mmax_index_list_1z]
		vzmax_list_1z = vz[mmax_index_list_1z]
		xmax_list_1z = x[mmax_index_list_1z]
		ymax_list_1z = y[mmax_index_list_1z]
		zmax_list_1z = z[mmax_index_list_1z]
		rvirmax_list_1z = rvir[mmax_index_list_1z]
		bmax_list_1z = b[mmax_index_list_1z]
		cmax_list_1z = c[mmax_index_list_1z]
		DAT = np.column_stack((IDnew_list_1z, mmax_list_1z, vxmax_list_1z, vymax_list_1z, vzmax_list_1z, xmax_list_1z, ymax_list_1z, zmax_list_1z, rvirmax_list_1z, bmax_list_1z, cmax_list_1z))
		np.savetxt('z_=%.4f_mmax.dat' % z_redshift, DAT, delimiter=" ")
		
		os.chdir("..")
		IDold = IDnew_list_1z
	
	np.savetxt('mmax_prog_deadendnum.dat', deadend_num_list, delimiter=" ")
	
def hist(var, pmax, pmin, binnum):
	binlist = np.arange(pmin, pmax + (pmax - pmin)/float(binnum), (pmax - pmin)/float(binnum))
#	xlocs = np.arange(pmin + .5 * (pmax - pmin)/binnum, pmax + .5 * (pmax - pmin)/binnum, (pmax - pmin)/binnum)
	histx, bin_edgesx = np.histogram(var, bins=binlist, normed = True)
	return bin_edgesx, histx

def percentile(var_list, i):
	var_list = np.sort(var_list)
	list_perc = [.68,.95, .997]
	N = len(var_list)
	upper_perc = .5 + list_perc[i-1] / 2e0
	lower_perc = .5 - list_perc[i-1] / 2e0
	return var_list[round(lower_perc * N)], var_list[round(upper_perc * N)]

def plotread_host(filelist0, l, variabletag):
#	FILE: AID_sat_1z, Acosangle_sat_1z, Av_relhost_cosangle_sat_1z, Av_relcent_cosangle_sat_1z, Apv_relhost_cosangle_sat_1z, Apv_relcent_cosangle_sat_1z, Adist_sat_1z, Am_sat_1z, Av_relhost_sat_1z

	if variabletag == 'progenitor_num':
		prog_num = loadtxt(filelist0_prog[l])
		pmax = 200
		pmin = 1
		var_list = prog_num
	if variabletag == 'logmass':
		mmax_list_1z = loadtxt(filelist0[l], usecols = (1, ))
#		print(np.log10(max(mmax_list_1z)), np.log10(min(mmax_list_1z)))
		pmax = 13.
		pmin = 8.
		var_list = np.log10(mmax_list_1z)
		
	elif variabletag == 'vel':
		vxmax_list_1z, vymax_list_1z, vzmax_list_1z = loadtxt(filelist0[l], usecols = (2,3,4), unpack = 1)
		threeDvel_1z = threeDvel(vxmax_list_1z, vymax_list_1z, vzmax_list_1z)
		pmax = 1500.
		pmin = 0.
		var_list = threeDvel_1z
		
	else:
		xmax_list_1z, ymax_list_1z, zmax_list_1z = loadtxt(filelist0[l], usecols = (5,6,7), unpack = 1)
		vxmax_list_1z, vymax_list_1z, vzmax_list_1z = loadtxt(filelist0[l], usecols = (2,3,4), unpack = 1)
		relvcosangle = []
		vcosangle = []
		ABdist = []
		A_mass = []
		B_mass = []
		for i in range(len(xmax_list_1z) - np.mod(len(xmax_list_1z), 2)):
			if i % 2 ==0:
				Apos = np.array([xmax_list_1z[i], ymax_list_1z[i], zmax_list_1z[i]])
				Bpos = np.array([xmax_list_1z[i+1], ymax_list_1z[i+1], zmax_list_1z[i+1]])		
				lineAB = coordsubstract(Bpos, Apos)/coorddistance(Apos, Bpos)
				lineBA = coordsubstract(Apos, Bpos)/coorddistance(Apos, Bpos)
				Avel = np.array([vxmax_list_1z[i], vymax_list_1z[i], vzmax_list_1z[i]])
				Bvel = np.array([vxmax_list_1z[i+1], vymax_list_1z[i+1], vzmax_list_1z[i+1]])
				Avel_norm = Avel/np.linalg.norm(Avel)
				Bvel_norm = Bvel/np.linalg.norm(Bvel)
				
				ABavg_vel = (Avel + Bvel)/2e0
				Avel_min_ABavg_norm = (Avel - ABavg_vel)/np.linalg.norm((Avel - ABavg_vel))
				Bvel_min_ABavg_norm = (Bvel - ABavg_vel)/np.linalg.norm((Bvel - ABavg_vel))
				
				relvcosangle.append(np.dot(Avel_min_ABavg_norm, lineAB))
				relvcosangle.append(np.dot(Bvel_min_ABavg_norm, lineBA))
				vcosangle.append(np.dot(Avel_norm, lineAB))
				vcosangle.append(np.dot(Bvel_norm, lineBA))
				ABdist.append(coorddistance(Apos, Bpos))
				
		if variabletag == 'dsep':
			pmax = 10e0
			pmin = 0e0
			var_list = ABdist
			
		if variabletag == 'cosangle_vel_to_AB':
			pmax = 1e0
			pmin = -1e0
			var_list = vcosangle

		if variabletag == 'cosangle_relvel_to_AB':
			pmax = 1e0
			pmin = -1e0
			var_list = relvcosangle
			
		if variabletag == 'b':
			bmax_list_1z = loadtxt(filelist0[l], usecols = (9,), unpack = 1)
			pmax = 1e0
			pmin = 0e0
			var_list = bmax_list_1z
	
	return pmax, pmin, var_list
	
def merger_analysis_plot2(variabletag):
#	mmax_list_1z, vxmax_list_1z, vymax_list_1z, vzmax_list_1z, xmax_list_1z, ymax_list_1z, zmax_list_1z, rvirmax_list_1z, bmax_list_1z, cmax_list_1z
	plot_list_mat = []
	plot_list = []
	plot_redshift_list = [z_redshift_list[0]]
	
	counter = 1
	var_list0 = np.array([])
	
	for l in range(len(filelist0)):
#		pmax, pmin, var_list = plotread_host(filelist0, l, variabletag)
#		print(z_redshift_list[l], 10e0**var_list[1000])
		
		if z_redshift_list[l] < z_redshift_list_linear[counter]:
			pmax, pmin, var_list = plotread_host(filelist0, l, variabletag)
#			print(z_redshift_list[l], max(var_list), min(var_list))
#			sys.exit()
			var_list0 = np.append(var_list0, var_list)
		else:
			var_list0.flatten()
			xlocs, histx = hist(var_list0, pmax, pmin, 50.)
			
			plot_list_mat.append(histx[::-1])
			plot_list.append([amax(var_list), amin(var_list), median(var_list), std(var_list)])	

			counter += 1
			var_list0 = np.array([])
			if l != len(filelist0)-1:
				plot_redshift_list.append(z_redshift_list[l])
				pmax, pmin, var_list = plotread_host(filelist0,l,variabletag)
				var_list0 = np.append(var_list0, var_list)
	xmin = 0
	xmax = z_redshift_list[-1]

	fig, ax = plt.subplots(1, 1, figsize=(6,6))
	ax.set_xlabel("z",fontsize=18)
	ax.set_xlim(xmin, xmax)
	ax.set_ylabel(ylabel[varname_list.index(variabletag)], fontsize=18)	

	plot_list_mat = np.transpose(plot_list_mat)
	im = plt.imshow(plot_list_mat, extent=[xmin, xmax, pmin, pmax], cmap=mpl.cm.rainbow, aspect = 'auto', interpolation='nearest')
	t = np.arange(amin(plot_list_mat.flatten()), amax(plot_list_mat.flatten()), (amax(plot_list_mat.flatten()) - amin(plot_list_mat.flatten())) / 5e0)
	im.set_clim(vmin=amin(plot_list_mat.flatten()), vmax=amax(plot_list_mat.flatten()))
	fig.colorbar(im, ticks=t, format='$%.2f$',cax = fig.add_axes([0.91, 0.15, 0.02, 0.7]))
	fig.savefig('traceMaxProg_host_' + variabletag + '.jpg', dpi=500)
	plt.close()

def merger_analysis_plot(variabletag):
#	mmax_list_1z, vxmax_list_1z, vymax_list_1z, vzmax_list_1z, xmax_list_1z, ymax_list_1z, zmax_list_1z, rvirmax_list_1z, bmax_list_1z, cmax_list_1z
	plot_list_mat = []
	plot_list = []
	lower_sigma_bounds = []
	upper_sigma_bounds = []
	mhost_max_list = []
	mhost_min_list = []
	z_redshift_list_sav = []
	z_mat_list = np.arange(0, )
	for l in range(27):
		z_redshift_list_sav.append(z_redshift_list[l])
	
	for l in range(27):
		if variabletag == 'progenitor_num':
			prog_num = loadtxt(filelist0_prog[l])
			pmax = 200
			pmin = 1
			var_list = prog_num
		if variabletag == 'logmass':
			mmax_list_1z = loadtxt(filelist0[l], usecols = (1, ))
			pmax = 13.
			pmin = 8.
			var_list = np.log10(mmax_list_1z)
			mhost_max_list.append(max(mmax_list_1z))
			mhost_min_list.append(min(mmax_list_1z))
			
			
		elif variabletag == 'vel':
			vxmax_list_1z, vymax_list_1z, vzmax_list_1z = loadtxt(filelist0[l], usecols = (2,3,4), unpack = 1)
			threeDvel_1z = threeDvel(vxmax_list_1z, vymax_list_1z, vzmax_list_1z)
			pmax = 1500.
			pmin = 0.
			var_list = threeDvel_1z
			
		else:
			xmax_list_1z, ymax_list_1z, zmax_list_1z = loadtxt(filelist0[l], usecols = (5,6,7), unpack = 1)
			vxmax_list_1z, vymax_list_1z, vzmax_list_1z = loadtxt(filelist0[l], usecols = (2,3,4), unpack = 1)
			relvcosangle = []
			vcosangle = []
			ABdist = []
			for i in range(len(xmax_list_1z) - np.mod(len(xmax_list_1z), 2)):
				if i % 2 ==0:
					Apos = np.array([xmax_list_1z[i], ymax_list_1z[i], zmax_list_1z[i]])
					Bpos = np.array([xmax_list_1z[i+1], ymax_list_1z[i+1], zmax_list_1z[i+1]])		
					lineAB = coordsubstract(Bpos, Apos)/coorddistance(Apos, Bpos)
					lineBA = coordsubstract(Apos, Bpos)/coorddistance(Apos, Bpos)
					Avel = np.array([vxmax_list_1z[i], vymax_list_1z[i], vzmax_list_1z[i]])
					Bvel = np.array([vxmax_list_1z[i+1], vymax_list_1z[i+1], vzmax_list_1z[i+1]])
					Avel_norm = Avel/np.linalg.norm(Avel)
					Bvel_norm = Bvel/np.linalg.norm(Bvel)
					
					ABavg_vel = (Avel + Bvel)/2e0
					Avel_min_ABavg_norm = (Avel - ABavg_vel)/np.linalg.norm((Avel - ABavg_vel))
					Bvel_min_ABavg_norm = (Bvel - ABavg_vel)/np.linalg.norm((Bvel - ABavg_vel))
					
					relvcosangle.append(np.dot(Avel_min_ABavg_norm, lineAB))
					relvcosangle.append(np.dot(Bvel_min_ABavg_norm, lineBA))
					vcosangle.append(np.dot(Avel_norm, lineAB))
					vcosangle.append(np.dot(Bvel_norm, lineBA))
					ABdist.append(coorddistance(Apos, Bpos))
					
			if variabletag == 'dsep':
				pmax = 10e0
				pmin = 0e0
				var_list = ABdist
				
			if variabletag == 'cosangle_vel_to_AB':
				pmax = 1e0
				pmin = -1e0
				var_list = vcosangle

			if variabletag == 'cosangle_relvel_to_AB':
				pmax = 1e0
				pmin = -1e0
				var_list = relvcosangle
				
			if variabletag == 'b':
				bmax_list_1z = loadtxt(filelist0[l], usecols = (9,), unpack = 1)
				pmax = 1e0
				pmin = 0e0
				var_list = bmax_list_1z
	
		xlocs, histx = hist(var_list, pmax, pmin, 100)
		plot_list_mat.append(histx[::-1])
		plot_list.append([np.amax(var_list), np.amin(var_list), np.median(var_list), np.std(var_list)])
		lower_1sigma, upper_1sigma = percentile(var_list, 1)
		lower_2sigma, upper_2sigma = percentile(var_list, 2)
		lower_3sigma, upper_3sigma = percentile(var_list, 3)
		lower_sigma_bounds.append([lower_1sigma, lower_2sigma, lower_3sigma])
		upper_sigma_bounds.append([upper_1sigma, upper_2sigma, upper_3sigma])
	
	if variabletag == 'logmass':
		DAT = np.column_stack((z_redshift_list_sav, mhost_max_list, mhost_min_list))
		np.savetxt('maxmin_hostmass.dat', DAT, delimiter=" ")

	xmin = 0
	xmax = z_redshift_list_sav[-1]
#	print(xmax)
	
	lower_sigma_bounds = np.transpose(lower_sigma_bounds)
	upper_sigma_bounds = np.transpose(upper_sigma_bounds)
		
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13,6))
	ax1.set_xlabel("z",fontsize=18)	
	ax1.set_xlim(xmin, xmax)
	ax1.set_ylim(pmin, pmax)
	ax1.set_ylabel(ylabel[varname_list.index(variabletag)], fontsize=18)	
	ax1.grid(True,linestyle= '-',which='major',color= '0.75')
	ax1.grid(True,linestyle= '-',which='minor',color= '0.75')
	ax1.grid(True, which='both')
	ax1.minorticks_on()
	plot_list = np.transpose(plot_list)
	
	ax1.plot(z_redshift_list_sav, plot_list[0])
	ax1.plot(z_redshift_list_sav, plot_list[1])
	ax1.fill_between(z_redshift_list_sav, plot_list[1], plot_list[0], facecolor='blue', alpha=0.2) # between max and min
	ax1.plot(z_redshift_list_sav, plot_list[2])
	ax1.fill_between(z_redshift_list_sav, lower_sigma_bounds[0], upper_sigma_bounds[0], facecolor='yellow', alpha=0.5) # one sigma
	ax1.fill_between(z_redshift_list_sav, lower_sigma_bounds[1], upper_sigma_bounds[1], facecolor='green', alpha=0.5) # two sigma
	ax1.fill_between(z_redshift_list_sav, lower_sigma_bounds[2], upper_sigma_bounds[2], facecolor='red', alpha=0.5) # three sigma
	
	ax2.set_ylabel(ylabel[varname_list.index(variabletag)], fontsize=18)
	
	plot_list_mat = np.transpose(plot_list_mat)
	print(plot_list_mat)
	plot_list_mat = np.array(plot_list_mat)
	im = plt.imshow(plot_list_mat, extent=[xmin, xmax, pmin, pmax], cmap=mpl.cm.rainbow, aspect = 'auto', interpolation='nearest')
#	if  variabletag == 'vel':
#		print(max(plot_list_mat.flatten()))
	t = np.arange(amin(plot_list_mat.flatten()), amax(plot_list_mat.flatten()), (amax(plot_list_mat.flatten())-amin(plot_list_mat.flatten()))/5e0)
	im.set_clim(vmin=amin(plot_list_mat.flatten()), vmax=amax(plot_list_mat.flatten()))
	
	ax2.set_xlabel("z",fontsize=18)
#	ax2.set_xticklabels(z_redshift_list, rotation='vertical', minor = 1)

#	ax2.xaxis.set_major_locator(majorLocator)
#	ax2.xaxis.set_major_formatter(majorFormatter)
	if variabletag == 'vel':
		fig.colorbar(im, ticks=t, format='$%.3f$',cax = fig.add_axes([0.91, 0.15, 0.02, 0.7]))
	else:
		fig.colorbar(im, ticks=t, format='$%.2f$',cax = fig.add_axes([0.91, 0.15, 0.02, 0.7]))
	
	fig.savefig('traceMaxProg_host_' + variabletag + 'old.jpg', dpi=500)
	plt.close()
	
if __name__ == "__main__":
##	 gathering all the catalog files into a list, ordered by snapshots from current day snapshots to beginning of the universe, 		# (scale a, scale = 1 is the same as redshift z = 0)
	Rockstarfilelist0 = []
	Rockstar_files = os.listdir(Rockstar_hlist_address)
	pwd0 = os.getcwd()

	for i in range(len(Rockstar_files)):
		filename = Rockstar_files[i]
		if filename.endswith("list"):
			Rockstarfilelist0.append(filename)
	Rockstarfilelist0.sort(key = natural_keys)
	Rockstarfilelist = Rockstarfilelist0[:][::-1]
	
	datmainfoldr = "/LGB_merger"
	os.chdir(pwd0  + datmainfoldr)

	merger_analysis()
	sys.exit()
	
	filelist0 = []
	filelist0_prog = []
	for file_1z in os.listdir('./'):
		if file_1z.endswith("_mmax.dat"):
			filelist0.append(file_1z)
#		if file_1z.endswith("_mmax_prognum.dat"):
#			filelist0_prog.append(file_1z)
	filelist0.sort(key = natural_keys)
	filelist0_prog.sort(key = natural_keys)

	z_redshift_list = []
	for l in range(len(filelist0)):
		z_redshift = float(filelist0[l][3:8])
		z_redshift_list.append(z_redshift)
		
	z_redshift_list_linear = np.arange(0, z_redshift_list[-1], .3)
	
#	ylabel = ['Number of Progenitor', '$\log(M_{host}$)', '$d_{sep}$', 'velocity', 'cosangle_vel_to_AB', 'cosangle_relvel_to_AB', 'b/a']
	ylabel = ['$\log(M_{host}$)', '$d_{sep}$', 'cosangle_vel_to_AB', 'cosangle_relvel_to_AB', 'b/a']
#	varname_list = ['progenitor_num', 'logmass', 'dsep', 'vel', 'cosangle_vel_to_AB', 'cosangle_relvel_to_AB', 'b']
	varname_list = ['logmass', 'dsep', 'cosangle_vel_to_AB', 'cosangle_relvel_to_AB', 'b']
#	for vartag in varname_list:
#	for vartag in ['progenitor_numb']:
	for vartag in ['logmass']:
#		merger_analysis_plot(vartag)
		merger_analysis_plot2(vartag)
	
