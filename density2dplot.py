import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter, MaxNLocator, AutoMinorLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from numpy import random, linspace, pi
import numpy as np
from matplotlib import rcParams
import sys
from histplot import *
from matplotlib.colors import Normalize
from helpers import *
import matplotlib as mpl

def densityplot_DsatVSDsep(x, y, xname, yname, plotpath, myvmin, myvmax):

	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.set_size_inches(5.5, 5.0, forward=True)
	axTemperature = fig.add_subplot(111)

	axTemperature.set_xlabel(xname,fontsize=14)
	axTemperature.set_ylabel(yname,fontsize=14)
	
	xmin, xmax, ymin, ymax = min(x), max(x), min(y), max(y)

	xlims = [xmin,xmax]
	ylims = [ymin,ymax]
	
	nxbins, nybins, nbins = 81, 81, 81
	
	xbins = linspace(start = xmin, stop = xmax, num = nxbins)
	ybins = linspace(start = ymin, stop = ymax, num = nybins)
	xcenter = (xbins[0:-1] + xbins[1:])/2.0
	ycenter = (ybins[0:-1] + ybins[1:])/2.0
	H, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins), normed = True)
		
	axTemperature.yaxis.set_major_locator(MaxNLocator(5))
	axTemperature.xaxis.set_major_locator(MaxNLocator(5))	
	axTemperature.xaxis.set_minor_locator(AutoMinorLocator(10))
	axTemperature.yaxis.set_minor_locator(AutoMinorLocator(8))
	
	plt.tick_params(axis='xy', colors='green')
	for tick in axTemperature.get_xticklines():
	    tick.set_color('black')
	for minortick in axTemperature.xaxis.get_minorticklines():
		minortick.set_color('black')	    
	for tick in axTemperature.get_yticklines():
	    tick.set_color('white')
	for minortick in axTemperature.yaxis.get_minorticklines():
		minortick.set_color('white')    

	axTemperature.text(1.25, .7, "normalized probability", color='black', fontsize=13, fontweight='bold', rotation=270, transform=axTemperature.transAxes)	

	# Plot the probability
	im = axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', aspect='auto', cmap= "Greys")
	im.set_clim(vmin=myvmin, vmax=myvmax)
	t = np.arange(0,myvmax+.3,.3)
	fig.colorbar(im, ticks=t, format='$%.2f$')
	
	axTemperature.set_xlabel(xname,fontsize=14)
	axTemperature.set_ylabel(yname,fontsize=14)

#	axTemperature.tick_params(axis='both', labelsize=10)
	#Set up the plot limits
	axTemperature.set_xlim(xlims)
	axTemperature.set_ylim(ylims)
	
	plt.tight_layout()
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500)
	plt.close()

def densityplot(x, y, xname, yname, plotpath, myvmin, myvmax):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.set_size_inches(4.5, 4.0, forward=True)
	axTemperature = fig.add_subplot(111)

	axTemperature.set_xlabel(xname,fontsize=14)
	axTemperature.set_ylabel(yname,fontsize=14)
	
	xmin, xmax, ymin, ymax = min(x), max(x), min(y), max(y)
	xlims = [xmin,xmax]
	ylims = [ymin,ymax]
	
	nxbins, nybins, nbins = 40, 40, 40
	
	xbins = linspace(start = xmin, stop = xmax, num = nxbins)
	ybins = linspace(start = ymin, stop = ymax, num = nybins)
	xcenter = (xbins[0:-1] + xbins[1:])/2.0
	ycenter = (ybins[0:-1] + ybins[1:])/2.0
	H, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins), normed = True)
		
	axTemperature.yaxis.set_major_locator(MaxNLocator(5))
	axTemperature.xaxis.set_major_locator(MaxNLocator(5))	
	axTemperature.xaxis.set_minor_locator(AutoMinorLocator(10))
	axTemperature.yaxis.set_minor_locator(AutoMinorLocator(8))
	
	plt.tick_params(axis='xy', colors='green')
	for tick in axTemperature.get_xticklines():
	    tick.set_color('white')
	for minortick in axTemperature.xaxis.get_minorticklines():
		minortick.set_color('white')	    
	for tick in axTemperature.get_yticklines():
	    tick.set_color('white')
	for minortick in axTemperature.yaxis.get_minorticklines():
		minortick.set_color('white')    

	# Plot the probability
	im = axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', aspect='auto')
	im.set_clim(vmin=myvmin, vmax=myvmax)
	t = np.arange(0,myvmax+.3,.3)
	fig.colorbar(im, ticks=t, format='$%.2f$')
	
	axTemperature.set_xlabel(xname,fontsize=14)
	axTemperature.set_ylabel(yname,fontsize=14)

#	axTemperature.tick_params(axis='both', labelsize=10)
	#Set up the plot limits
	axTemperature.set_xlim(xlims)
	axTemperature.set_ylim(ylims)
	
	plt.tight_layout()
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500)
	plt.close()
		
def densityplot_wsides(x, y, xname, yname, plotpath, limx):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.set_size_inches(4.5, 4.0, forward=True)
	gs = gridspec.GridSpec(3, 3)
	gs.update(wspace=0, hspace=0)
	axHistx = plt.subplot(gs[0, :-1])
	axTemperature = plt.subplot(gs[1:, :-1])
	axHisty = plt.subplot(gs[1:, 2])

	# Find the min/max of the data
	xmin = limx[0]
	xmax = limx[1]
	ymin = 0e0
	ymax = 5e-1

	xlims = [xmin,xmax]
	ylims = [ymin,ymax]
	
	nxbins = 40
	nybins = 40
	nbins = 40
	
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	xlocs = np.arange(xmin + .5*(xmax-xmin)/nbins, xmax + .5*(xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ylocs = np.arange(ymin + .5*(ymax-ymin)/nbins, ymax + .5*(ymax-ymin)/nbins, (ymax-ymin)/nbins)	
		
	xleft = np.delete(xbins, [40])
	yleft = np.delete(ybins, [40])
	axHistx.set_xlim(xmin, xmax)
	axHisty.set_ylim(ymin, ymax)
	axHisty.set_xlim(0, 4)
	axHistx.set_ylim(.3, 1.)
	#Plot the histograms
	histx, bin_edgesx = np.histogram(x, bins=xbins, normed = True)
	histy, bin_edgesy = np.histogram(y, bins=ybins, normed = True)

	axHistx.plot(xlocs, histx, color='blue')
	axHisty.plot(histy, ylocs, color='red')
	axHistx.fill_between(xlocs, histx-0.0030813452641, histx+0.0030813452641, facecolor='red', alpha=0.2)
		
	#Make the tickmarks pretty
	ticklabels = axHistx.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(10)
	ticklabels = axHisty.get_xticklabels()
	for label in ticklabels:
		label.set_fontsize(10)

	# Remove the inner axes numbers of the histograms
	axHistx.xaxis.set_major_formatter(NullFormatter())
	axHistx.xaxis.set_ticks_position('none') 
	axHisty.yaxis.set_major_formatter(NullFormatter())
	axHisty.yaxis.set_ticks_position('none') 
	axHistx.yaxis.set_major_locator(MaxNLocator(4))		
	axHisty.xaxis.set_major_locator(MaxNLocator(4))
	
#	xbins = linspace(start = xmin, stop = xmax, num = nxbins)
#	ybins = linspace(start = ymin, stop = ymax, num = nybins)

	H, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins), normed = True)

	# Plot the probability
	im = axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', aspect='auto')
	t = np.arange(0,H.max(),H.max()/5)
	fig.colorbar(im, ticks=t, format='$%.2f$')
	im.set_clim(vmin=0, vmax=3)
	axTemperature.yaxis.set_major_locator(MaxNLocator(5, prune='upper'))
	axTemperature.xaxis.set_major_locator(MaxNLocator(5, prune='upper'))
	axTemperature.xaxis.set_minor_locator(AutoMinorLocator(10))
	axTemperature.yaxis.set_minor_locator(AutoMinorLocator(8))
		
	axTemperature.set_xlabel(xname,fontsize=14)
	axTemperature.set_ylabel(yname,fontsize=14)

	axTemperature.tick_params(axis='both', labelsize=10)
	
	#Set up the plot limits
	axTemperature.set_xlim(xlims)
	axTemperature.set_ylim(ylims)
	
	# Save to a File	
	plt.tight_layout()

	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500)
	plt.close()

def densityplot_w1sides(x, y, xname, yname, plotpath, limx):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.set_size_inches(4.5, 4.0, forward=True)
	gs = gridspec.GridSpec(6,6)
	gs.update(wspace=0, hspace=0)
	axHistx = plt.subplot(gs[:1, :-2])
	axTemperature = plt.subplot(gs[1:, :-1])

	# Find the min/max of the data
	xmin = limx[0]
	xmax = limx[1]
	ymin = 0e0
	ymax = 5e-1

	xlims = [xmin,xmax]
	ylims = [ymin,ymax]
	
	nxbins = 40
	nybins = 40
	nbins = 40
	
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	xlocs = np.arange(xmin + .5*(xmax-xmin)/nbins, xmax + .5*(xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ylocs = np.arange(ymin + .5*(ymax-ymin)/nbins, ymax + .5*(ymax-ymin)/nbins, (ymax-ymin)/nbins)	
		
	xleft = np.delete(xbins, [40])
	yleft = np.delete(ybins, [40])
	axHistx.set_xlim(xmin, xmax)
	axHistx.set_ylim(.3, 1.)
	#Plot the histograms
	histx, bin_edgesx = np.histogram(x, bins=xbins, normed = True)
	histy, bin_edgesy = np.histogram(y, bins=ybins, normed = True)
	axHistx.plot(xlocs, histx, color='blue')
	axHistx.fill_between(xlocs, histx - 0.0030813452641, histx + 0.0030813452641, facecolor='red', alpha=0.2)
	
	#Make the tickmarks pretty
	ticklabels = axHistx.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(10)

	# Remove the inner axes numbers of the histograms
	axHistx.xaxis.set_major_formatter(NullFormatter())
	axHistx.xaxis.set_ticks_position('none') 
	axHistx.yaxis.set_major_locator(MaxNLocator(4))		

	H, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins), normed = True)

	# Plot the probability
	im = axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', aspect='auto')
#	t = np.arange(0,H.max(),H.max()/5)
	t = np.arange(0, 3.3, .3)
	fig.colorbar(im, ticks=t, format='$%.2f$')
	axTemperature.yaxis.set_major_locator(MaxNLocator(5, prune='upper'))
	axTemperature.xaxis.set_major_locator(MaxNLocator(5, prune='upper'))
	axTemperature.xaxis.set_minor_locator(AutoMinorLocator(10))
	axTemperature.yaxis.set_minor_locator(AutoMinorLocator(8))
		
	axTemperature.set_xlabel(xname,fontsize=14)
	axTemperature.set_ylabel(yname,fontsize=14)

	axTemperature.tick_params(axis='both', labelsize=10)
	im.set_clim(vmin=0, vmax=3)
	#Set up the plot limits
	axTemperature.set_xlim(xlims)
	axTemperature.set_ylim(ylims)

	plt.tight_layout()
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500)
	plt.close()
	
def densityplot_2panels_signal(x1, y1, x2, y2, xname, yname, plotpath, limx, limy, lim_hist, hist2dlim_max, binno):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
	fig = plt.figure()
	fig.set_size_inches(6.5, 4.5, forward=True)
	
	plt.subplots_adjust(hspace=0.0)
	gsmom = gridspec.GridSpec(1, 2, width_ratios = [.95, .95], hspace = .0, wspace = .0) 
	#make nested gridspecs
	gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = gsmom[0], height_ratios = [1, 2])
	gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = gsmom[1], height_ratios = [1, 2])
	
#	gs = gridspec.GridSpec(6,12)
	
#	gs.update(left=0.1, right=0.9, wspace=0.0)
	axHistx1 = plt.subplot(gs1[0])
	axTemperature1 = plt.subplot(gs1[1])
	axHistx2 = plt.subplot(gs2[0])
	axTemperature2 = plt.subplot(gs2[1])
#	axHistx3 = plt.subplot(gs2[1])
#	axTemperature3 = plt.subplot(gs2[3])
	# Find the min/max of the data
	
	xmin = limx[0]
	xmax = limx[1]
	ymin = limy[0]
	ymax = limy[1]
	hist_min = lim_hist[0]
	hist_max = lim_hist[1]
	
	axHistx1.set_ylabel("$p(\cos \\theta)$",fontsize=16)
	xlims = [xmin,xmax]
	ylims = [ymin,ymax]
	
	nxbins = binno
	nybins = binno
	nbins = binno
#	
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	xlocs = np.arange(xmin + .5*(xmax-xmin)/nbins, xmax + .5*(xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ylocs = np.arange(ymin + .5*(ymax-ymin)/nbins, ymax + .5*(ymax-ymin)/nbins, (ymax-ymin)/nbins)	
	
	xleft = np.delete(xbins, [nbins])
	yleft = np.delete(ybins, [nbins])
	ax_x_list = [axHistx1, axHistx2]
	ax_T_list = [axTemperature1, axTemperature2]	
	xlist = [x1,x2]
	ylist = [y1,y2]
	
	Qminus10 = []
	Qplus90 = []
	rls = []
	histx_12 = []
	for i in range(len(xlist)):
		ax_x_list[i].set_xlim(xmin, xmax)
		ax_x_list[i].set_ylim(hist_min, hist_max)
		#Plot the histograms
		histx, bin_edgesx = np.histogram(xlist[i], bins=xbins)
		histx = histx/sum(histx) * binno/2
		histx_12.append(histx)
		
#		Qminus10_0 = sum(histx[:int(binno*0.5)])
#		Qplus90_0 = sum(histx[int(binno*0.5):])

		Qminus10_0 = sum(histx[:int(binno*0.1)])
		Qplus90_0 = sum(histx[int(binno*0.9):])
				
		Qminus10.append(Qminus10_0)
		Qplus90.append(Qplus90_0)
		
		meanAB, sigmaAB = markovN(len(xlist[i]))
		
		rls.append((meanAB, sigmaAB))
				
		ax_x_list[i].plot(xlocs, histx, color='blue')
	
		#Make the tickmarks pretty

		ticklabels = ax_x_list[i].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(18)
		ticklabels = ax_x_list[i].get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(18)	
		ticklabels = ax_T_list[i].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(18)
		ticklabels = ax_T_list[i].get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(18)	
					
		# Remove the inner axes numbers of the histograms
		ax_x_list[i].xaxis.set_major_formatter(NullFormatter())
		ax_x_list[i].yaxis.set_major_locator(MaxNLocator(4, prune='lower'))		

		H, xedges,yedges = np.histogram2d(ylist[i],xlist[i],bins=(ybins,xbins), normed = True)
		# Plot the probability
		im = ax_T_list[i].imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', aspect='auto', cmap = 'jet', vmin=0, vmax=hist2dlim_max)
		ax_T_list[i].set_xlabel(xname,fontsize=18)

		if i ==0:
			ax_T_list[i].set_ylabel(yname,fontsize=18)
		ax_T_list[i].xaxis.set_major_locator(MaxNLocator(5, prune='upper'))
#		ax_T_list[i].yaxis.set_major_locator(MaxNLocator(5, prune='upper'))
		
		ax_T_list[i].xaxis.set_minor_locator(AutoMinorLocator(10))
		ax_T_list[i].yaxis.set_minor_locator(AutoMinorLocator(8))
		if i ==0:
			ax_T_list[i].tick_params(axis='both')
		else:
			ax_T_list[i].yaxis.set_ticklabels([])
			ax_x_list[i].yaxis.set_ticklabels([])
		ax_T_list[i].set_xlim(xlims)
		ax_T_list[i].set_ylim(ylims)
		
	ax_x_list[0].set_title('(a) Signal', fontsize = 23)
	ax_x_list[1].set_title('(b) Control', fontsize = 23)
#	ax_x_list[2].set_title('(c) Signal/Control', fontsize = 23)

	i = 0
#	histx_ratio = (histx_12[i] / histx_12[i+1])/sum(histx_12[i] / histx_12[i+1]) * binno/2e0
#	
#	Qminus10.append(sum(histx_ratio[:int(binno*0.1)]))
#	Qplus90.append(sum(histx_ratio[int(binno*0.9):]))

##	Qminus10.append(sum(histx_ratio[:int(binno*0.5)]))
##	Qplus90.append(sum(histx_ratio[int(binno*0.5):]))
#	
#	meanAB, sigmaAB = markovN(len(xlist[1]))
#	
#	rls.append((meanAB, sigmaAB))
	
#	ax_x_list[2].set_xlim(xmin, xmax)
#	ax_x_list[2].set_ylim(hist_min, hist_max)
#	ax_x_list[2].yaxis.set_major_locator(MaxNLocator(4, prune='both'))
#	ax_x_list[2].plot(xlocs, np.array(histx_ratio), color='blue')
#	
#	ticklabels = ax_x_list[2].get_yticklabels()
#	for label in ticklabels:
#		label.set_fontsize(9)

#	# Remove the inner axes numbers of the histograms
#	ax_x_list[2].xaxis.set_major_formatter(NullFormatter())
#	ax_x_list[2].yaxis.set_ticklabels([])

#	H, xedges, yedges = np.histogram2d(ylist[i], xlist[i], bins=(ybins, xbins))
#	H2, xedges2, yedges2 = np.histogram2d(ylist[i + 1], xlist[i + 1], bins=(ybins, xbins))
#	
#	H_ratio = []
#	H_sum = 0
#	for j in range(len(H)):
#		H_ratio.append(np.array([]))
#		for k in range(len(H[j])):
##			if H2[j][k] != 0:
#			H_ratio[j] = np.append( H_ratio[j], H[j][k]/H2[j][k])
#				H_sum += H[j][k] - H2[j][k]
#			else:
#				H_ratio[j] = np.append(H_ratio[j], 0)

#	for j in range(len(H)):
#		for k in range(len(H[j])):
#			H_ratio[j][k] = binno ** 2e0 * H_ratio[j][k]/H_sum
	
#	ax_T_list[2].yaxis.set_ticklabels([])
	
	for j in range(1):
		print(Qplus90[j],Qminus10[j], rls[j])
		ax_x_list[j].text(-.3, .6, "$q = %.3f , %.0f \sigma$" % (Qplus90[j]/Qminus10[j], abs((Qplus90[j]/Qminus10[j] - rls[j][0])/rls[j][1])), color='red', fontsize=16, fontweight='bold')
		
#	ax_T_list[2].set_xlabel(xname,fontsize=14)

#	ax_T_list[2].xaxis.set_major_locator(MaxNLocator(5))
#	ax_T_list[2].xaxis.set_minor_locator(AutoMinorLocator(10))
#	ax_T_list[2].yaxis.set_minor_locator(AutoMinorLocator(8))
#	ax_T_list[2].yaxis.set_ticklabels([])
#	ax_T_list[2].set_xlim(xlims)
#	ax_T_list[2].set_ylim(ylims)
	
#	im = ax_T_list[2].imshow(H_ratio, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', cmap = 'jet', aspect='auto', norm=Normalize(vmin=0.6, vmax=1.7))
	ax_T_list[-1].text(1.2, 1, "normalized probability", color='black', fontsize=16, fontweight='bold', rotation=270, transform=ax_T_list[-1].transAxes)		
	t = np.arange(0,hist2dlim_max,0.3)
	cb = plt.colorbar(im, ticks=t, format='$%.1f$',cax = fig.add_axes([0.91, 0.15, 0.02, 0.7]))
	cb.ax.tick_params(labelsize=13)
	im.set_clim(vmin=0, vmax=hist2dlim_max)
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500,bbox_inches='tight')
	plt.close()
	
def densityplot_3panels(x1, y1, x2, y2, x3, y3, xname, yname, plotpath, plottitles, limx, limy, lim_hist, hist2dlim_max, binno):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.set_size_inches(10, 4.5, forward=True)
	
	plt.subplots_adjust(hspace=0.0)
	gsmom = gridspec.GridSpec(1, 2, width_ratios = [.95, 2], hspace = .0, wspace = .0) 
	#make nested gridspecs
	gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = gsmom[0], height_ratios = [1, 2])
	gs2 = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec = gsmom[1], wspace = .0, hspace = .0, height_ratios = [1, 2], width_ratios = [1,1])
	
#	gs = gridspec.GridSpec(6,12)
	
#	gs.update(left=0.1, right=0.9, wspace=0.0)
	axHistx1 = plt.subplot(gs1[0])
	axTemperature1 = plt.subplot(gs1[1])
	axHistx2 = plt.subplot(gs2[0])
	axTemperature2 = plt.subplot(gs2[2])
	axHistx3 = plt.subplot(gs2[1])
	axTemperature3 = plt.subplot(gs2[3])
	# Find the min/max of the data
	
	xmin = limx[0]
	xmax = limx[1]
	ymin = limy[0]
	ymax = limy[1]
	hist_min = lim_hist[0]
	hist_max = lim_hist[1]
	
	axHistx1.set_ylabel("$p(\cos \\theta)$",fontsize=18)
	xlims = [xmin,xmax]
	ylims = [ymin,ymax]
	
	nxbins = binno
	nybins = binno
	nbins = binno
#	
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	xlocs = np.arange(xmin + .5*(xmax-xmin)/nbins, xmax + .5*(xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ylocs = np.arange(ymin + .5*(ymax-ymin)/nbins, ymax + .5*(ymax-ymin)/nbins, (ymax-ymin)/nbins)	
		
	xleft = np.delete(xbins, [nbins])
	yleft = np.delete(ybins, [nbins])
	ax_x_list = [axHistx1, axHistx2, axHistx3]
	ax_T_list = [axTemperature1, axTemperature2, axTemperature3]	
	xlist = [x1,x2,x3]
	ylist = [y1,y2,y3]
	
	for i in range(len(xlist)):
		ax_x_list[i].set_xlim(xmin, xmax)
		ax_x_list[i].set_ylim(hist_min, hist_max)
		#Plot the histograms
		histx, bin_edgesx = np.histogram(xlist[i], bins=xbins, normed = True)
		
		indpos = np.where(xlist[i] > 0.8)
		indneg = np.where(xlist[i] < -0.8)
		Qplus90 = len(indpos[0])
		Qminus10 = len(indneg[0])
		
		meanAB, sigmaAB = markovN(len(xlist[i]))

		ax_x_list[i].plot(xlocs, histx, color='blue')
	
		#Make the tickmarks pretty
		ticklabels = ax_x_list[i].get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(18)		
		ticklabels = ax_x_list[i].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(18)

		ticklabels = ax_T_list[i].get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(18)		
		ticklabels = ax_T_list[i].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(18)
			
		# Remove the inner axes numbers of the histograms
		ax_x_list[i].xaxis.set_major_formatter(NullFormatter())
		ax_x_list[i].yaxis.set_major_locator(MaxNLocator(4, prune='lower'))		

		H, xedges,yedges = np.histogram2d(ylist[i],xlist[i],bins=(ybins,xbins), normed = True)
		# Plot the probability
		im = ax_T_list[i].imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', cmap = 'jet', origin='lower', aspect='auto', vmin=0, vmax=hist2dlim_max)
		ax_T_list[i].set_xlabel(xname,fontsize=14)
		if i ==0:
			ax_T_list[i].set_ylabel(yname,fontsize=14)
		ax_T_list[i].xaxis.set_major_locator(MaxNLocator(5, prune='upper'))
#		ax_T_list[i].yaxis.set_major_locator(MaxNLocator(5, prune='upper'))
		ticklabels = ax_T_list[i].get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(13)
		ticklabels = ax_T_list[i].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(13)
								
		ax_T_list[i].xaxis.set_minor_locator(AutoMinorLocator(10))
		ax_T_list[i].yaxis.set_minor_locator(AutoMinorLocator(8))
		if i ==0:
			ax_T_list[i].tick_params(axis='both')
		else:
			ax_T_list[i].yaxis.set_ticklabels([])
			ax_x_list[i].yaxis.set_ticklabels([])
		ax_T_list[i].set_xlim(xlims)
		ax_T_list[i].set_ylim(ylims)
		ax_x_list[i].text(-.1, .6, "$q = %.2f, %.f \sigma$" % (Qplus90/Qminus10, (Qplus90/Qminus10 - meanAB)/sigmaAB), color='red', fontsize=13, fontweight='bold')	

	ax_x_list[0].set_title(plottitles[0], fontsize = 23)
	ax_x_list[1].set_title(plottitles[1], fontsize = 23)
	ax_x_list[2].set_title(plottitles[2], fontsize = 23)	
	
	ax_T_list[-1].text(1.3, 1, "normalized probability", color='black', fontsize=14, fontweight='bold', rotation=270, transform=ax_T_list[-1].transAxes)		
	t = np.arange(0,hist2dlim_max,0.3)
	cbar = plt.colorbar(im, ticks=t, format='$%.2f$', cax = fig.add_axes([0.91, 0.15, 0.02, 0.7]))
	cbar.ax.tick_params(labelsize=14) 
	im.set_clim(vmin=0, vmax=hist2dlim_max)
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500,bbox_inches='tight')
	plt.close()

		
def densityplot_2panels(x1, y1, x2, y2, xname, yname, plotpath, limx, limy, lim_hist, hist2dlim_max, binno):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.set_size_inches(6, 4.8, forward=True)
	gs = gridspec.GridSpec(3,6)
	gs.update(left=0.1, right=0.9, wspace=0.1)
	axHistx1 = plt.subplot(gs[0:1, :3])
	axTemperature1 = plt.subplot(gs[1:, :3])
	axHistx2 = plt.subplot(gs[0:1,3:6])
	axTemperature2 = plt.subplot(gs[1:, 3:6])
	# Find the min/max of the data
	
	axHistx1.set_ylabel("$p(\cos \\theta)$",fontsize=14)
	xmin = limx[0]
	xmax = limx[1]
	ymin = limy[0]
	ymax = limy[1]
	hist_min = lim_hist[0]
	hist_max = lim_hist[1]
	
	xlims = [xmin,xmax]
	ylims = [ymin,ymax]
	
	nxbins = binno
	nybins = binno
	nbins = binno
#	
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	xlocs = np.arange(xmin + .5*(xmax-xmin)/nbins, xmax + .5*(xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ylocs = np.arange(ymin + .5*(ymax-ymin)/nbins, ymax + .5*(ymax-ymin)/nbins, (ymax-ymin)/nbins)	
		
	xleft = np.delete(xbins, [nbins])
	yleft = np.delete(ybins, [nbins])
	ax_x_list = [axHistx1, axHistx2]
	ax_T_list = [axTemperature1, axTemperature2]	
	xlist = [x1,x2]
	ylist = [y1,y2]

	for i in range(len(ax_x_list)):
		ax_x_list[i].set_xlim(xmin, xmax)
		ax_x_list[i].set_ylim(hist_min, hist_max)
		#Plot the histograms
		histx, bin_edgesx = np.histogram(xlist[i], bins=xbins)
		histy, bin_edgesy = np.histogram(ylist[i], bins=ybins)
		histx = histx/sum(histx) * binno/2
		histy = histy/sum(histy) * binno/2
		
		indpos = np.where(xlist[i] > 0.8)
		indneg = np.where(xlist[i] < -0.8)
		Qplus90 = len(indpos[0])
		Qminus10 = len(indneg[0])
		
		meanAB, sigmaAB = markovN(len(xlist[i]))
				
		ax_x_list[i].plot(xlocs, histx, color='blue')
	
		# Remove the inner axes numbers of the histograms
		ax_x_list[i].xaxis.set_major_formatter(NullFormatter())
		ax_x_list[i].yaxis.set_major_locator(MaxNLocator(4))		

		ticklabels = ax_x_list[i].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(15)
		ticklabels = ax_x_list[i].get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(15)	
		ticklabels = ax_T_list[i].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(15)
		ticklabels = ax_T_list[i].get_xticklabels()
		for label in ticklabels:
			label.set_fontsize(15)
			
		H, xedges,yedges = np.histogram2d(ylist[i],xlist[i],bins=(ybins,xbins), normed = True)
		# Plot the probability
		im = ax_T_list[i].imshow(H, extent=[xmin,xmax,ymin,ymax], cmap = 'jet', interpolation='nearest', origin='lower', aspect='auto', vmin = 0, vmax = hist2dlim_max)
		ax_T_list[i].set_xlabel(xname,fontsize=14)
		
		if i ==0:
			ax_T_list[i].set_ylabel(yname,fontsize=14)
#		t = np.arange(0,H.max(),H.max()/5)
	
		ax_T_list[i].yaxis.set_major_locator(MaxNLocator(5, prune='upper'))
		if i ==0:
			ax_T_list[i].xaxis.set_major_locator(MaxNLocator(5, prune='upper'))
		else:
			ax_T_list[i].xaxis.set_major_locator(MaxNLocator(5))		
		ax_T_list[i].xaxis.set_minor_locator(AutoMinorLocator(10))
		ax_T_list[i].yaxis.set_minor_locator(AutoMinorLocator(8))
		if i ==0:
			ax_T_list[i].tick_params(axis='both')
		else:
			ax_T_list[i].yaxis.set_ticklabels([])
			ax_x_list[i].yaxis.set_ticklabels([])
		ax_T_list[i].set_xlim(xlims)
		ax_T_list[i].set_ylim(ylims)

		ax_x_list[i].text(-.1, .6, "$q = %.3f, %.0f \sigma$" % (Qplus90/Qminus10, (Qplus90/Qminus10 - meanAB)/sigmaAB), color='red', fontsize=13, fontweight='bold')
	
	t = np.arange(0,hist2dlim_max+.3,.3)
	cb = plt.colorbar(im, ticks=t, format='$%.1f$',cax = fig.add_axes([0.91, 0.15, 0.02, 0.7]))
	cb.ax.tick_params(labelsize=15)
	
	ax_T_list[-1].text(1.25, 1, "normalized probability", color='black', fontsize=14, fontweight='bold', rotation=270, transform=ax_T_list[-1].transAxes)
	im.set_clim(vmin=0, vmax=hist2dlim_max)
#	plt.tight_layout()
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500,bbox_inches='tight')
	plt.close()
	
def densityratioplot_2panels(x1, y1, x2, y2, x3, y3, xname, yname, plotpath, limx, binno):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.set_size_inches(6, 4.5, forward=True)
	gs = gridspec.GridSpec(3,6)
	gs.update(left=0.1, right=0.9, wspace=0.1)
	axHistx1 = plt.subplot(gs[0:1, :3])
	axTemperature1 = plt.subplot(gs[1:, :3])
	axHistx2 = plt.subplot(gs[0:1,3:6])
	axTemperature2 = plt.subplot(gs[1:, 3:6])
	# Find the min/max of the data
	xmin = limx[0]
	xmax = limx[1]
	ymin = 0e0
	ymax = 5e-1
	axHistx1.set_ylabel("$p(\cos \\theta)$",fontsize=14)
	xlims = [xmin,xmax]
	ylims = [ymin,ymax]
	
	nxbins = binno
	nybins = binno
	nbins = binno
	
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	xlocs = np.arange(xmin + .5e0*(xmax-xmin)/nbins, xmax + .5e0*(xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ylocs = np.arange(ymin + .5e0*(ymax-ymin)/nbins, ymax + .5e0*(ymax-ymin)/nbins, (ymax-ymin)/nbins)	
		
	xleft = np.delete(xbins, [nbins])
	yleft = np.delete(ybins, [nbins])
	
	ax_x_list = [axHistx1, 0, axHistx2]
	ax_T_list = [axTemperature1, 0, axTemperature2]	
	
	xlist = [x1,x2,x1,x3]
	ylist = [y1,y2,y1,y3]

	histxlist = []
	for i in range(4):
		histx, bin_edgesx = np.histogram(xlist[i], bins=xbins, normed = True)
		histxlist.append(histx)
	
	for i in [0,2]:
		histx = histxlist[i]
		histx2 = histxlist[i+1]
		histx_ratio = []
		for j in range(len(histx)):
			histx_ratio.append(histx[j] / histx2[j])
		ax_x_list[i].set_xlim(xmin, xmax)
		ax_x_list[i].set_ylim(.5, 1.4)			
		ax_x_list[i].plot(xlocs, histx_ratio, color='blue')
		ticklabels = ax_x_list[i].get_yticklabels()
		for label in ticklabels:
			label.set_fontsize(10)

		# Remove the inner axes numbers of the histograms
		ax_x_list[i].xaxis.set_major_formatter(NullFormatter())
		ax_x_list[i].yaxis.set_major_locator(MaxNLocator(4))		

		# Plot the probability
		if i == 0:
			ax_T_list[i].set_ylabel(yname,fontsize=14)
#	t = np.arange(0,H.max(),H.max()/5)
	
		ax_T_list[i].yaxis.set_major_locator(MaxNLocator(5, prune='upper'))
		if i == 0:
			ax_T_list[i].xaxis.set_major_locator(MaxNLocator(5, prune='upper'))
		else:
			ax_T_list[i].xaxis.set_major_locator(MaxNLocator(5))	
		ax_T_list[i].xaxis.set_minor_locator(AutoMinorLocator(10))
		ax_T_list[i].yaxis.set_minor_locator(AutoMinorLocator(8))
		ax_T_list[i].set_xlabel(xname,fontsize=14)
		if i ==0:
			ax_T_list[i].tick_params(axis='both')
		else:
			ax_T_list[i].yaxis.set_ticklabels([])
			ax_x_list[i].yaxis.set_ticklabels([])
		ax_T_list[i].set_xlim(xlims)
		ax_T_list[i].set_ylim(ylims)
		
		H, xedges, yedges = np.histogram2d(ylist[i], xlist[i], bins=(ybins,xbins), normed = True)
		H2, xedges2, yedges2 = np.histogram2d(ylist[i+1], xlist[i+1], bins=(ybins,xbins), normed = True)
#		print (H[0])
		H_ratio = []
		for j in range(len(H)):
			H_ratio.append([])
			for k in range(len(H[j])):
				H_ratio[j].append(H[j][k]/ H2[j][k])
		im = ax_T_list[i].imshow(H_ratio, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', aspect='auto', norm=Normalize(vmin=0.6, vmax=1.7))
	
	t = np.arange(0.6,1.7,.1)
	plt.colorbar(im, ticks=t, format='$%.2f$',cax = fig.add_axes([0.91, 0.15, 0.02, 0.7]))
	im.set_clim(vmin=0.6, vmax=1.7)
	ax_T_list[-1].text(1.25, 1, "normalized probability", color='black', fontsize=14, fontweight='bold', rotation=270, transform=ax_T_list[-1].transAxes)
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500,bbox_inches='tight')
	plt.close()

def densityplot_8panels(xlistall, ylistall, xname, yname, featurename, plotpath, limx, binno):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
	fig = plt.figure()
	fig.set_size_inches(8, 13, forward=True)
	gs = gridspec.GridSpec(13,4)
	gs.update(left=0.1, right=0.9, hspace=0.01, wspace=0.01)
	axHistxlist = []
	axTemperaturelist = []
	for j in range(4):
		axHistxlist.append([])
		axTemperaturelist.append([])
		for i in range(2):
			axHistxlist[j].append(plt.subplot(gs[j * 3 : j * 3 + 1, i : i + 1]))
			axTemperaturelist[j].append(plt.subplot(gs[j * 3 + 1 : j * 3 + 3, i : i + 1]))
			
	xmin, xmax, ymin,ymax = -1e0, 1e0, 0e0, .5e0

	xlims, ylims = [xmin,xmax], [ymin,ymax]
	
	nxbins, nybins, nbins = binno, binno, binno
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	xlocs = np.arange(xmin + .5*(xmax-xmin)/nbins, xmax + .5*(xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ylocs = np.arange(ymin + .5*(ymax-ymin)/nbins, ymax + .5*(ymax-ymin)/nbins, (ymax-ymin)/nbins)	
		
	xleft = np.delete(xbins, [nbins])
	yleft = np.delete(ybins, [nbins])
	subfigname = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)']
	for j in range(4):
		xlist = xlistall[j]
		ylist = ylistall[j]
		axHistxlist[j][0].set_ylabel("$p(\cos \\theta)$",fontsize=14)
		for i in range(2):
			axHistxlist[j][i].set_xlim(xmin, xmax)
			axHistxlist[j][i].set_ylim(.4, .85)
			#Plot the histograms
			histx, bin_edgesx = np.histogram(xlist[i], bins=xbins, normed = True)
			axHistxlist[j][i].plot(xlocs, histx, color='blue')
			axHistxlist[j][i].text(0.42, 0.6, subfigname[i+j*2], color='black', fontsize=20, fontweight='bold', transform=axHistxlist[j][i].transAxes)
			#Make the tickmarks pretty
			ticklabels = axHistxlist[j][i].get_yticklabels()
			for label in ticklabels:
				label.set_fontsize(10)

			# Remove the inner axes numbers of the histograms
			axHistxlist[j][i].xaxis.set_major_formatter(NullFormatter())
			axHistxlist[j][i].yaxis.set_major_locator(MaxNLocator(4))		

			H, xedges,yedges = np.histogram2d(ylist[i],xlist[i],bins=(ybins,xbins), normed = True)
			# Plot the probability
			im = axTemperaturelist[j][i].imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', aspect='auto', vmin=0, vmax=3)
			if j == 3:
				axTemperaturelist[j][i].set_xlabel(xname,fontsize=14)
			if i == 0:
				axTemperaturelist[j][i].set_ylabel(yname,fontsize=14)
				if j != 2:
					axTemperaturelist[j][i].text(-.85, 0.5, featurename[j], color='red', fontsize=20, fontweight='bold', transform=axTemperaturelist[j][i].transAxes)
				else:
					axTemperaturelist[j][i].text(-.85, 0.5, featurename[j], color='red', fontsize=18, fontweight='bold', transform=axTemperaturelist[j][i].transAxes)
			axTemperaturelist[j][i].yaxis.set_major_locator(MaxNLocator(5, prune='upper'))
			if i != 1:
				axTemperaturelist[j][i].xaxis.set_major_locator(MaxNLocator(5, prune='upper'))
			else:
				axTemperaturelist[j][i].xaxis.set_major_locator(MaxNLocator(5))		
			axTemperaturelist[j][i].xaxis.set_minor_locator(AutoMinorLocator(10))
			axTemperaturelist[j][i].yaxis.set_minor_locator(AutoMinorLocator(8))
			if i ==0:
				axTemperaturelist[j][i].tick_params(axis='both')
			else:
				axTemperaturelist[j][i].yaxis.set_ticklabels([])
				axHistxlist[j][i].yaxis.set_ticklabels([])
			if j != 3:
				axTemperaturelist[j][i].xaxis.set_ticklabels([])
			if i != 0:
				axHistxlist[j][i].yaxis.set_ticklabels([])				
			axTemperaturelist[j][i].set_xlim(xlims)
			axTemperaturelist[j][i].set_ylim(ylims)

	t = np.arange(0,3.3,.3)
	plt.colorbar(im, ticks=t, format='$%.2f$',cax = fig.add_axes([0.61, 0.2, 0.02, 0.6]))
	axTemperaturelist[1][-1].text(1.15, .4, "normalized probability", color='black', fontsize=20, fontweight='bold', rotation=270, transform=axTemperaturelist[1][-1].transAxes)	
	im.set_clim(vmin=0, vmax=3)
#	plt.tight_layout(pad=1.0)
#	gs.tight_layout(fig, rect=[0.1, 0, .9, 1])
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500,bbox_inches='tight')
#	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500,bbox_inches='tight')	
	plt.close()

def densityplot_MERGER_SAT(xlistall, ylistall, xnames, ynames, redshiftnames, xlims, ylims, plotpath, binno,plottag):
	fig = plt.figure()
	fig.set_size_inches(14.5, 16, forward=True)
	gs = gridspec.GridSpec(5,4)
	gs.update(left=0.1, right=0.9, hspace=0.15, wspace=0.35)
	axTemperaturelist = []
	for j in range(len(xlistall)):
		axTemperaturelist.append([])
		for i in range(len(xnames)):
			axTemperaturelist[j].append(plt.subplot(gs[j:j+1, i:i+1]))
	
	nbins = binno
#	subfigname = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']
	for j in range(len(xlistall)):
		xlist = xlistall[j]
		ylist = ylistall[j]
		for i in range(len(xnames)):
			[xmin, xmax] = xlims[i]
			[ymin, ymax] = ylims[i]
			
			xbins = np.linspace(xmin, xmax, nbins+1)
			ybins = np.linspace(ymin, ymax, nbins+1)
			xlocs = np.arange(xmin + .5*(xmax-xmin)/nbins, xmax + .5*(xmax-xmin)/nbins, (xmax-xmin)/nbins)
			ylocs = np.arange(ymin + .5*(ymax-ymin)/nbins, ymax + .5*(ymax-ymin)/nbins, (ymax-ymin)/nbins)
			X, Y = np.meshgrid(xlocs, ylocs)
			
			H, yedges, xedges = np.histogram2d(ylist[i], xlist[i], bins=(ybins, xbins))
#			Hx, bin_edgesx = np.histogram(xlist[i], bins=nbins, normed = True)
#			Hy, bin_edgesy = np.histogram(ylist[i], bins=nbins, normed = True)
#			Hx = np.meshgrid(Hx, Hy)
			axTemperaturelist[j][i].contour(X, Y, H, colors='black')
#			print(xbins, ybins, H)
			
			# Plot the probability
			if plottag == 'hostA':
				vmax0 = 1700
			if plottag == 'hostB':
				vmax0 = 1200			
			im = axTemperaturelist[j][i].imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', aspect='auto', vmin=0, vmax=vmax0)
			if j == len(xlistall)-1:
				axTemperaturelist[j][i].set_xlabel(xnames[i], fontsize=25)
			axTemperaturelist[j][i].set_ylabel(ynames[i],fontsize=25)
			if i == 0:
				axTemperaturelist[j][i].text(-0.5, 0.75, redshiftnames[j], color='red', fontsize=15, fontweight='bold', transform=axTemperaturelist[j][i].transAxes)

			axTemperaturelist[j][i].xaxis.set_major_locator(MaxNLocator(4))		
			axTemperaturelist[j][i].xaxis.set_minor_locator(AutoMinorLocator(8))
			axTemperaturelist[j][i].yaxis.set_major_locator(MaxNLocator(5))	
			axTemperaturelist[j][i].yaxis.set_minor_locator(AutoMinorLocator(10))
#			if i ==0:
#				axTemperaturelist[j][i].tick_params(axis='both')
#			else:
#				axTemperaturelist[j][i].yaxis.set_ticklabels([])
			if j != len(xlistall)-1:
				axTemperaturelist[j][i].xaxis.set_ticklabels([])
			axTemperaturelist[j][i].set_xlim(xlims[i])
			axTemperaturelist[j][i].set_ylim(ylims[i])
			
	if plottag == 'hostA':
		t = np.arange(0, 2000, 300)
		im.set_clim(vmin=0, vmax=2000)
	if plottag == 'hostB':
		t = np.arange(0, 1200, 300)	
		im.set_clim(vmin=0, vmax=900)
	plt.colorbar(im, ticks=t, format='$%.0f$',cax = fig.add_axes([0.91, 0.2, 0.02, 0.6]))
	axTemperaturelist[1][-1].text(1.5, .4, "number per bin", color='black', fontsize=15, fontweight='bold', rotation=270, transform=axTemperaturelist[1][-1].transAxes)	
	
#	plt.tight_layout(pad=1.0)
#	gs.tight_layout(fig, rect=[0.1, 0, .9, 1])
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500, bbox_inches='tight')
	plt.close()

def densityplot_MERGER_SAT_NORM(xlistall, ylistall, xnames, ynames, redshiftnames, xlims, ylims, plotpath, nbins, plottag):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.set_size_inches(25.5, 6, forward=True)
	gs = gridspec.GridSpec(1,4)
	gs.update(left=0, right=1.0, hspace=0., wspace=0.)
	axTemperaturelist = []
	for j in range(len(xnames)):
		axTemperaturelist.append([])
		for i in range(len(xlistall)):
			axTemperaturelist[j].append(plt.subplot(gs[j:j+1, i:i+1]))
				
#	subfigname = ['(a)', '(b)', '(c)']
	for i in range(len(xlistall)):
		xlist = xlistall[i]
		ylist = ylistall[i]
		for j in range(len(xnames)):
			if i == 0:
				ylims = [0, 0.5]
			else:
				ylims = [0, 1.7]
			[xmin, xmax] = xlims[0]
			[ymin, ymax] = ylims

			xbins = np.linspace(xmin, xmax, nbins+1)
			ybins = np.linspace(ymin, ymax, nbins+1)
			xlocs = np.arange(xmin + .5*(xmax-xmin)/nbins, xmax + .5*(xmax-xmin)/nbins, (xmax-xmin)/nbins)
			ylocs = np.arange(ymin + .5*(ymax-ymin)/nbins, ymax + .5*(ymax-ymin)/nbins, (ymax-ymin)/nbins)
			X, Y = np.meshgrid(xlocs, ylocs)
			
			H, yedges, xedges = np.histogram2d(ylist[j], xlist[j], bins=(ybins, xbins), normed = True)
			axTemperaturelist[j][i].contour(X, Y, H, colors='black')
			
			myaspect = (xmax-xmin)/  (ymax-ymin)
			im = axTemperaturelist[j][i].imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', cmap = 'jet', aspect=myaspect, vmin = 0, vmax = 1.6)
			
			axTemperaturelist[j][i].set_ylabel(ynames[j], fontsize=25)
			axTemperaturelist[j][i].set_xlabel(xnames[j],fontsize=25)
			axTemperaturelist[j][i].text(0.35, 1.05, redshiftnames[i], color='red', fontsize=30, fontweight='bold', transform=axTemperaturelist[j][i].transAxes)

			axTemperaturelist[j][i].xaxis.set_major_locator(MaxNLocator(4))		
			axTemperaturelist[j][i].xaxis.set_minor_locator(AutoMinorLocator(8))
			axTemperaturelist[j][i].yaxis.set_major_locator(MaxNLocator(5))	
			axTemperaturelist[j][i].yaxis.set_minor_locator(AutoMinorLocator(10))
			
			ticklabels = axTemperaturelist[j][i].get_yticklabels()
			for label in ticklabels:
				label.set_fontsize(22)
			ticklabels = axTemperaturelist[j][i].get_xticklabels()
			for label in ticklabels:
				label.set_fontsize(22)
						
#			if i ==0:
#				axTemperaturelist[j][i].tick_params(axis='both')
#			else:
#				axTemperaturelist[j][i].yaxis.set_ticklabels([])
#			if j != len(xlistall)-1:
#				axTemperaturelist[j][i].xaxis.set_ticklabels([])
			axTemperaturelist[j][i].set_xlim(xlims[j])
			axTemperaturelist[j][i].set_ylim(ylims[j])
			if i != 0:
				axTemperaturelist[j][i].plot([-1,-1],[0,0.5], c = 'red', linewidth = 5)
				axTemperaturelist[j][i].plot([-1,1],[0.5,0.5], c = 'red', linewidth = 5)
				axTemperaturelist[j][i].plot([1,1],[0.5,0], c = 'red', linewidth = 5)
				axTemperaturelist[j][i].plot([1,-1],[0,0], c = 'red', linewidth = 5)
			
	axTemperaturelist[0][-1].text(1.25, .85, "normalized probability", color='black', fontsize=29, fontweight='bold', rotation=270, transform=axTemperaturelist[0][-1].transAxes)	
	cb = plt.colorbar(im, format='$%.1f$')
	cb.ax.tick_params(labelsize=25)
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500, bbox_inches='tight')
	plt.close()

def densityplot_MERGER_SAT_NORM2(xlistall, ylistall, xnames, ynames, redshiftnames, xlims, ylims, plotpath, nbins, plottag):
	fig = plt.figure()
	fig.set_size_inches(35, 10, forward=True)
	gs = gridspec.GridSpec(2,6)
	gs.update(left=0.1, right=0.8, hspace=0.1, wspace=0.1)
	axTemperaturelist = []
	for i in range(len(xlistall)):
		axTemperaturelist.append([])
		for j in range(len(xnames)):
			axTemperaturelist[i].append(plt.subplot(gs[j:j+1, i:i+1]))
			
	subfigname = ['(a) Flyby', '(b) Non-flyby']
	for i in range(len(xlistall)):
		xlist = xlistall[i]
		ylist = ylistall[i]
		for j in range(len(xnames)):
			[xmin, xmax] = xlims[j]
			[ymin, ymax] = ylims[j]
			
			xbins = np.linspace(xmin, xmax, nbins+1)
			ybins = np.linspace(ymin, ymax, nbins+1)
			xlocs = np.arange(xmin + .5 * (xmax-xmin)/nbins, xmax + .5 * (xmax-xmin)/nbins, (xmax - xmin)/nbins)
			ylocs = np.arange(ymin + .5 * (ymax-ymin)/nbins, ymax + .5 * (ymax-ymin)/nbins, (ymax - ymin)/nbins)
			X, Y = np.meshgrid(xlocs, ylocs)
			
			H, yedges, xedges = np.histogram2d(ylist[j], xlist[j], bins=(ybins, xbins), normed = True)
			axTemperaturelist[i][j].contour(X, Y, H, colors='black')
			
			myaspect = (xmax-xmin)/  (ymax-ymin)
			im = axTemperaturelist[i][j].imshow(H, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', aspect=myaspect)
			if i == 0:
				axTemperaturelist[i][j].set_ylabel(ynames[j], fontsize=30)
#				axTemperaturelist[i][j].text(0.45, -0.35, subfigname[i], color='black', fontsize=20, fontweight='bold', transform=axTemperaturelist[i][j].transAxes)
				
			axTemperaturelist[i][j].set_xlabel(xnames[j],fontsize=30)
			if j == 0:
				axTemperaturelist[i][j].set_title(redshiftnames[i], color='red', fontsize=25, fontweight='bold')
#			if j == 0:
#				axTemperaturelist[i][j].set_title(subfigname[i], fontsize=23)
			axTemperaturelist[i][j].xaxis.set_major_locator(MaxNLocator(4))		
			axTemperaturelist[i][j].xaxis.set_minor_locator(AutoMinorLocator(8))
			axTemperaturelist[i][j].yaxis.set_major_locator(MaxNLocator(5))	
			axTemperaturelist[i][j].yaxis.set_minor_locator(AutoMinorLocator(10))
			plt.colorbar(im, ax = axTemperaturelist[i][j], format='$%.0f$')
			
			if j != len(xlistall)-1:
				axTemperaturelist[i][j].xaxis.set_ticklabels([])
			axTemperaturelist[i][j].set_xlim(xlims[j])
			axTemperaturelist[i][j].set_ylim(ylims[j])
		
	axTemperaturelist[-1][0].text(1.25, .4, "normalized probability", color='black', fontsize=25, fontweight='bold', rotation=270, transform=axTemperaturelist[-1][0].transAxes)
	
	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500, bbox_inches='tight')
	plt.close()
							
def densityratioplot_w1sides(x, y, x2, y2, xname, yname, plotpath):
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.set_size_inches(4.5, 4.0, forward=True)
	gs = gridspec.GridSpec(6,6)
	gs.update(wspace=0, hspace=0)
	axHistx = plt.subplot(gs[:1, :-2])
	axTemperature = plt.subplot(gs[1:, :-1])
#	axHisty = plt.subplot(gs[1:, 2])

	# Find the min/max of the data
	xmin = -1e0
	xmax = 1e0
	ymin = 0e0
	ymax = 5e-1

	xlims = [xmin,xmax]
	ylims = [ymin,ymax]
	
	nxbins = 40
	nybins = 40
	nbins = 40
	
	xbins = np.arange(xmin, xmax + (xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ybins = np.arange(ymin, ymax + (ymax-ymin)/nbins, (ymax-ymin)/nbins)
	xlocs = np.arange(xmin + .5*(xmax-xmin)/nbins, xmax + .5*(xmax-xmin)/nbins, (xmax-xmin)/nbins)
	ylocs = np.arange(ymin + .5*(ymax-ymin)/nbins, ymax + .5*(ymax-ymin)/nbins, (ymax-ymin)/nbins)	
	xleft = np.delete(xbins, [40])
	yleft = np.delete(ybins, [40])
	axHistx.set_xlim(xmin, xmax)
#	axHisty.set_ylim(ymin, ymax)
	
	#Plot the histograms
	histx, bin_edgesx = np.histogram(x, bins=xbins, normed = True)
	histy, bin_edgesy = np.histogram(y, bins=ybins, normed = True)
	histx2, bin_edgesx2 = np.histogram(x2, bins=xbins, normed = True)
	histy2, bin_edgesy2 = np.histogram(y2, bins=ybins, normed = True)
	
	histx_ratio = []
	histy_ratio = []
	for i in range(len(histx)):
		histx_ratio.append(histx2[i] / histx[i])
		histy_ratio.append(histy2[i] / histy[i])
	axHistx.set_ylim(0, 1.4)
#	axHisty.set_xlim(0, 1.2)		
#	axHisty.set_ylim(0, 0.5)
	axHistx.set_xlim(-1, 1)
			
	widthx = (xmax-xmin)/nbins
	widthy = (ymax-ymin)/nbins
	axHistx.plot(xlocs, histx_ratio, color='blue')
#	axHisty.plot(histy_ratio,ylocs, color='red')
	
	#Make the tickmarks pretty
	ticklabels = axHistx.get_yticklabels()
	for label in ticklabels:
		label.set_fontsize(10)
#	ticklabels = axHisty.get_xticklabels()
#	for label in ticklabels:
#		label.set_fontsize(10)

	# Remove the inner axes numbers of the histograms
	axHistx.xaxis.set_major_formatter(NullFormatter())
	axHistx.xaxis.set_ticks_position('none') 
#	axHisty.yaxis.set_major_formatter(NullFormatter())
#	axHisty.yaxis.set_ticks_position('none') 
	axHistx.yaxis.set_major_locator(MaxNLocator(3))		
#	axHisty.xaxis.set_major_locator(MaxNLocator(3))
	
	H, xedges,yedges = np.histogram2d(y,x,bins=(ybins,xbins), normed = True)
	H2, xedges2,yedges2 = np.histogram2d(y2,x2,bins=(ybins,xbins), normed = True)
	
	H_ratio = []
	for i in range(len(H)):
		H_ratio.append([])
		for j in range(len(H[i])):
			H_ratio[i].append(H2[i][j]/ H[i][j])
#	print len(H_ratio), len(H_ratio[0])
	# Plot the probability
	im = axTemperature.imshow(H_ratio, extent=[xmin,xmax,ymin,ymax], interpolation='nearest', origin='lower', aspect='auto', norm=Normalize(vmin=0.5, vmax=1.4))
#	axcolor = fig.add_axes(np.arange(0, 1.5, 0.1))
	t = np.arange(0,1.4,0.2)
	fig.colorbar(im, ticks=t, format='$%.2f$')
	
	axTemperature.yaxis.set_major_locator(MaxNLocator(5, prune='upper'))
	axTemperature.xaxis.set_major_locator(MaxNLocator(5, prune='upper'))
	axTemperature.xaxis.set_minor_locator(AutoMinorLocator(10))
	axTemperature.yaxis.set_minor_locator(AutoMinorLocator(8))
	
	axTemperature.set_xlabel(xname,fontsize=14)
	axTemperature.set_ylabel(yname,fontsize=14)
	
#	im.set_clim(vmin=0, vmax=3)
	axTemperature.tick_params(axis='both', labelsize=10)
	
	#Set up the plot limits
	axTemperature.set_xlim(xlims)
	axTemperature.set_ylim(ylims)
	
	# Save to a File	
	plt.tight_layout()

	plt.savefig(plotpath, format = 'pdf', transparent=True, dpi=500,bbox_inches='tight')
	plt.close()		

