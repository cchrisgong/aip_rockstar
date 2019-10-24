
'''
This code studies the merger tree history of the LGB of today z = 0
'''

from numpy import loadtxt, histogram
import numpy as np
import sys
from constant import *
import os
from shutil import copyfile
from helpers import *
import pandas as pd

def merger():
	# gathering all the catalog files into a list, ordered by snapshots from current day snapshots to beginning of the universe, 		# (scale a, scale = 1 is the same as redshift z = 0)
	Rockstarfilelist0 = []
	Rockstar_files = os.listdir(Rockstar_hlist_address)
	pwd0 = os.getcwd()

	for i in range(len(Rockstar_files)):
		filename = Rockstar_files[i]
		if filename.endswith("list"):
			Rockstarfilelist0.append(filename)
	Rockstarfilelist0.sort(key = natural_keys)
	
	Rockstarfilelist = Rockstarfilelist0[:][::-1]
	
	z0_filename = Rockstarfilelist[0] #z = 0 catalog file
	snapnum = z0_filename[6:12]
	
	datmainfoldr = "/LGB_merger"
	try:
		os.makedirs(datmainfoldr)
	except OSError:
		pass
			
	datfoldr = "/snapshot_for_z=%.4f" % 0.0
	
	datfoldrz0 = "/snapshot" + snapnum
	
	satfoldrz0 = pwd0 + datmainfoldr + datfoldr + "/smallhalos_highres" #folder of z = 0 sat files

	sat_files0 = os.listdir(satfoldrz0)
	sat_files0.sort(key = natural_keys)
	
	z0_LGB_path = pwd0 + datfoldrz0 + "/" + high_res_trueLG_pop #z = 0 LGB file
	

	os.chdir(pwd0  + "/" + datmainfoldr)

	LGBsat_path_old = satfoldrz0
	outputfile = "output_merger.dat" #output file for the entire merger calculation

	processornum = int(sys.argv[1])
#	processornum = 1
#	chunksize = 10
	chunksize = 215
#	for l in range(1, 2): #for each redshift outside z = 0
	for l in range(1, len(Rockstarfilelist)):
		filename_high_res = Rockstarfilelist[l]
		
		a = float(filename_high_res[6:12]) #scale
		z = (1.- a)/a #redshift
		
		# copy hlist file to home address
		filepath_catalog = Rockstar_hlist_address + filename_high_res
		newfilepath = os.getcwd() + "/snapshot_for_z=%.4f.list" % z
		if os.path.isfile(newfilepath): #if file exists
			pass
		else:
			copyfile(filepath_catalog, newfilepath) 
		
		# loading parent halo ID (from the last red shift) 
#		newID = loadtxt(newfilepath, skiprows = headerlength, usecols = (IDcol,), dtype='i8')#loading IDs of binaries at z= 0
		newdescID = loadtxt(newfilepath, skiprows = headerlength, usecols = (descIDcol,), dtype='i8')
		# for each redshift there is a folder to store data in
		datfoldr = "snapshot_for_z=%.4f" % z
		try:
			os.makedirs(datfoldr)
		except OSError:
			pass
		
		os.chdir(os.getcwd() + "/" + datfoldr)
		
		# merger tree satellites
		LGBsat_path_new = "smallhalos_highres_merg"

		f = open(newfilepath,"r") #read data
		lines = f.readlines()
		f.close()		

		try:
			os.makedirs(LGBsat_path_new)
		except OSError:
			pass

		for i in range((processornum - 1) * chunksize, min(processornum * chunksize, len( sat_files0 ))):
			old_satfile = LGBsat_path_old + "/" + sat_files0[i]
			old_satID = pd.read_csv(old_satfile, skiprows = 0, usecols = (IDcol,), delim_whitespace = True)
			old_satID = old_satID.values
			
			index = np.where( np.in1d( newdescID, old_satID ) )
			f_newsat = open(LGBsat_path_new + "/" + sat_files0[i], "w")
			f_newsat.write(lines[0])
			
#			if i == 2:
#				print("##########",newdescID[index[0]], old_satID) # checked		
	
			for j in range(len(index[0])):
				f_newsat.write(lines[index[0][j] + headerlength]) #checked
			f_newsat.close()
		
		pwd1 = os.getcwd()
		LGBsat_path_old = pwd1 + "/" + LGBsat_path_new

		os.chdir("..")
#		os.remove(newfilepath)
		
if __name__ == "__main__":
	merger()
	
