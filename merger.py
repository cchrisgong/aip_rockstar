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
	print(z0_filename)
	snapnum = z0_filename[6:12]
	datfoldrz0 = "/snapshot" + snapnum 
	satfoldrz0 = pwd0 + datfoldrz0 + "/smallhalos_highres" #folder of z = 0 sat files
	z0_LGB_path = pwd0 + datfoldrz0 + "/" + high_res_trueLG_pop #z = 0 LGB file
	
	IDold = loadtxt(z0_LGB_path, skiprows = 1, usecols = (IDcol,), dtype='i8')#loading IDs of binaries at z= 0
	
	datmainfoldr = "LGB_merger"
	try:
		os.makedirs(datmainfoldr)
	except OSError:
		pass
	os.chdir(pwd0  + "/" + datmainfoldr)
				
	LGBsat_path_old = satfoldrz0
	outputfile = "output_merger.dat" #output file for the entire merger calculation
	
	for l in range(1, len(Rockstarfilelist)): #for each redshift outside z = 0
		filename_high_res = Rockstarfilelist[l]
		print(filename_high_res)
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
		newID = pd.read_csv(newfilepath, skiprows = headerlength, sep = ' ', usecols = (IDcol,), dtype='i8')
		newdescID = pd.read_csv(newfilepath, skiprows = headerlength, sep = ' ', usecols = (descIDcol,), dtype='i8')		
		newID = newID.values
		newdescID = newdescID.values
#		newID, newdescID  = loadtxt(newfilepath, skiprows = 1, usecols = (IDcol, descIDcol), unpack=True, dtype='i8')
		# output file for the entire merger calculation
		f_o = open(outputfile, 'w')
		f_o.write("new ID unpacked for "+ str(l) + "/" + str( len( Rockstarfilelist ) ) + "\n")
		f_o.close()
		
		# for each redshift there is a folder to store data in
		datfoldr = "snapshot_for_z=%.4f" % z
		try:
			os.makedirs(datfoldr)
		except OSError:
			pass
		
		os.chdir(os.getcwd() + "/" + datfoldr)
		
		# merger tree local group binaries
		LGB_path_new = os.getcwd() + "/" + "LGB_z=%.4f.dat" % z
		f2 = open(LGB_path_new, "w")
#		
		f = open(newfilepath,"r") #read data
		lines = f.readlines()
		f.close()
#		
		f2.write(lines[0]) #write header
				
		index = np.where( np.in1d( newdescID, IDold ) )
#		print (newdescID[0], newdescID[headerlength], lines[headerlength], lines[2 * headerlength])
		for j in range(len(index[0])):
#			print(newdescID[index[0][j]], IDold[j])
#			print(lines[index[0][j] + headerlength + 1])
			f2.write(lines[index[0][j] + headerlength + 1])
		f2.close()

		IDold = newID[index[0]]
		
		#output file for the entire merger calculation
		f_o = open(outputfile, 'w')
		f_o.write("LGB done "+ str(l) + "/" + str(len(Rockstarfilelist)) + "\n")
		f_o.close()

		os.chdir("..")
#		os.remove(newfilepath)
		
		
if __name__ == "__main__":
	merger()
	
