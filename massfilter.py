'''
Massfilter.py shorten the data file snap_135_new.z0.000.AHF_halos by first elliminating the small mass halos. Mass range is 1e12 to 2e12 though it can also be slightly bigger to include more, such as 5e11 to 2e12. 

Massfilter has another function of picking out the super massive halos beyond this range. This is to further eliminate local groups that otherwise qualify for not being around any equally massive ones from 5e11 to 2e12, but are within range of say a third mass weighing 1e13. 
'''
from numpy import loadtxt
import numpy as np
import sys
from constant import *

def massfilter(LGmassmin, LGmassmax, OMmassmin, SMALLmassmax, filename, filenameLG_filtered, filenameOM_filtered, masscol, outputfile):
	# Read mass data from .AHF_halos file
	hrm = loadtxt(filename, skiprows = 1, usecols = (masscol,), unpack=True)	
	
	counterLG = 0
	counterOM = 0
	f2 = open(filenameLG_filtered,"w")
	f3 = open(filenameOM_filtered,"w")

	f = open(filename,"r") #read header
	line = f.readline()
	f.close()
	f2.write(line)# skip header
	f3.write(line)

	f = open(filename,"r") #read header
	lines = f.readlines()
	f.close()
	
	LGindex = np.where( ( hrm > LGmassmin ) & ( hrm < LGmassmax ) )
	OMindex = np.where( hrm > LGmassmax )
	for i in LGindex[0]:
		f2.write(lines[i+headerlength])
		counterLG += 1 #count the number of objects
		
	for i in OMindex[0]:
		f3.write(lines[i+headerlength])
		counterOM += 1
	f = open(outputfile, "a")
	f.write("ESMD simulation file: " + filename + "\n")
	f.write("number of objects with mass range (M > " + str(LGmassmin) + " and M < " + str(LGmassmax) + "): " + str(counterLG) + "\n")
	f.write("filtered LG mass halo entries are stored under: " + filenameLG_filtered + "\n")
	f.close()
	print ("ESMD simulation file: " + filename)
	print ("number of objects with mass range (M > " + str(LGmassmin) + " and M < " + str(LGmassmax) + "): " + str(counterLG))
	print ("filtered LG mass halo entries are stored under: " + filenameLG_filtered)
	print ("number of objects with mass range (M > " + str(OMmassmin) + ")" + str(counterOM))
	print ("filtered LG mass halo entries are stored under: " + filenameOM_filtered)
	f2.close()
	f3.close()
