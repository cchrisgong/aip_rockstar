'''
running on Rockstar higher resolution sample
'''

#import matplotlib.pyplot as plt
import sys
import os
from shutil import copyfile
from numpy import loadtxt
import numpy as np
from constant import *
from helpers import *
import pandas as pd

from massfilter import *
from massfunction import *
#from Lgfinder import *
from Lgfinderchecker import *
from checknosubhalo import *
from histplot import *
#from subhaloFinderIndexed import *
from nonindexed_smallhaloFinder import *
#from smallhalosNonindexedXYZPlot_rot import *
from smallhalosAngleDep import *
#from sieveLGgap import *

if __name__ == "__main__":
	Rockstarfilelist0 = []
	Rockstar_files = os.listdir(Rockstar_hlist_address)
	
	for i in range(len(Rockstar_files)):
		filename = Rockstar_files[i]
		if filename.endswith("list"):
			Rockstarfilelist0.append(filename)
	Rockstarfilelist0.sort(key = natural_keys)
	
	Rockstarfilelist0 = Rockstarfilelist0[30:][::-1]

	z_redshift_list = []
	Rockstarfilelist = []
	for i in range(len(Rockstarfilelist0)):
		if i % 4 == 0:
			Rockstarfilelist.append(Rockstarfilelist0[i])
			z_redshift_list.append(z_redshift[i])
	
	print(z_redshift)
	processornum = 1
	chunksize = 1
		
	for l in range((processornum - 1) * chunksize, processornum * chunksize):
		print(l)
		filename_high_res = Rockstarfilelist[l]
		
		snapnum = filename_high_res[6 : 12]
		datfoldr = "snapshot" + snapnum
		
		try:
			os.makedirs(datfoldr)
		except OSError:
			pass
		
		filepath_catalog = Rockstar_hlist_address + filename_high_res
		newfilepath = "./" + datfoldr + "/" + filename_high_res
		if os.path.isfile(newfilepath): #if file exists
			pass
		else:		
			copyfile(filepath_catalog, newfilepath)
		os.chdir("./" + datfoldr)
		outputfile = "output" + str(snapnum) + ".dat"

		LGmassmax = LGmassmax_z0
		LGmassmin = LGmassmin_z0
				
		f = open(outputfile, "w")
		f.write(filename_high_res + '\n' + "Mmin = " + str(LGmassmin) + '\n' + "z = " + str(z_redshift_list[l]) + '\n' + "Mmax = " + str(LGmassmax) + '\n')
		f.close()

		##### MASSFILTER ##########
		
		massfilter(LGmassmin, LGmassmax, OMmassmin, SMALLmassmax, filename_high_res, filename_high_res_filtered_lg_wsub, filename_high_res_filtered_overmassive, masscol, outputfile) #2 is the column for the mass
		print ("########################### HR MASS FILTER COMPLETED ########################\n")

		num_nosubs, num_justsubs = checknosubhalo(filename_high_res_filtered_lg_wsub, filename_high_res_filtered_lg_nosub, filename_high_res_filtered_lg_justsub) #throw away binaries components which are subhalos from the true LG halos binaries selection
		f = open(outputfile, "a")
		f.write("number of local group haloes without subs:" + str(num_nosubs) + '\n')
		f.write("number of local group haloes that are subs:" + str(num_justsubs) + '\n')
		f.close()

#		#LG mass range haloes info
		hrallLGID, hrallLGm, hrallLGx, hrallLGy, hrallLGz, hrallLGRvir = col_reader(filename_high_res_filtered_lg_nosub, IDcol, masscol, xcol, ycol, zcol, Rvircol) #no subs

		hrsubLGm, hrsubLGx, hrsubLGy, hrsubLGz = col_reader0(filename_high_res_filtered_lg_justsub, masscol, xcol, ycol, zcol) #only subs
		
		#overmassive halo info
		hromx, hromy, hromz = xyz_reader2(filename_high_res_filtered_overmassive, xcol, ycol, zcol)
				
#		####### LG SEARCH #########

#		#HIGH RES N^2 SEARCH
		LGbruteforceFinder(hrallLGm, hrallLGx, hrallLGy, hrallLGz, hromx, hromy, hromz, hrsubLGm, hrsubLGx, hrsubLGy, hrsubLGz, filename_high_res_bf, outputfile)	
		print ("############# HIGH RES LGB N^2 BRUTE FORCE SEARCH COMPLETED ##########\n")
		####### LG SEARCH #########
#		
#		
		
