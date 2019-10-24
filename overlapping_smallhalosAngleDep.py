from numpy import random
from scipy import stats
import os
from numpy import loadtxt, arctan2, sqrt, cos, fabs, sqrt, linspace
import matplotlib.pyplot as plt
import numpy as np
from helpers import *
from matplotlib.ticker import NullFormatter, MaxNLocator
import sys
from constant import *
#from overlapping import *
#from density2dplot import *
#from smallhalosAngleDep import *

'''
plotting the cos of angle distribution of the satellite galaxies to the host galaxies for overlapping signals
'''


def cos_satellite_angle(mA, mB, hostApos, hostBpos, Sepdis, SMALL):
	[SMALLm, SMALLx, SMALLy, SMALLz] = SMALL
	
	cosA = []
	cosB = []	
	
	#using half distance
	rA, rB = Sepdis/2e0, Sepdis/2e0
	
	for i in range(len(SMALLm)):
		Sx, Sy, Sz = SMALLx[i], SMALLy[i], SMALLz[i]
		pos = np.array([Sx, Sy, Sz])
		mass = SMALLm[i]
		satDistA = coorddistance(pos, hostApos)
		satDistB = coorddistance(pos, hostBpos)
		if satDistA < rA: #rA from host A
			cosangle = lawofcos(pos, hostApos, hostBpos)
			cosA.append([cosangle, mass, satDistA, Sx, Sy, Sz])

		if satDistB < rB: #rA from host A
			cosangle = lawofcos(pos, hostBpos, hostApos)		
			cosB.append([cosangle, mass, satDistB, Sx, Sy, Sz])
	return cosA, cosB
	
def cos_satellite_angle_OL(mA, mB, hostApos, hostBpos, Sepdis, SMALL):
	[SMALLx, SMALLy, SMALLz, SMALLm] = SMALL
	
	cosA = []
	cosB = []	
	
	#using half distance
	rA, rB = Sepdis/2e0, Sepdis/2e0
	
	for i in range(len(SMALLx)):
		Sx, Sy, Sz = SMALLx[i], SMALLy[i], SMALLz[i]
		Sm = SMALLm[i]
		pos = np.array([Sx, Sy, Sz])
		satDistA = coorddistance(pos, hostApos)
		satDistB = coorddistance(pos, hostBpos)
		if satDistA < rA: #rA from host A
			cosangle = lawofcos(pos, hostApos, hostBpos)
			cosA.append([cosangle, satDistA, Sm])

		if satDistB < rB: #rA from host A
			cosangle = lawofcos(pos, hostBpos, hostApos)		
			cosB.append([cosangle, satDistB, Sm])
	return cosA, cosB
	
def cosdistn_fake(LGID, LGm, LGx, LGy, LGz, LGRvir, LGb, LGc, sepDistlist, cosine_path, folder_control):
	#plotting histogram of the cosine angle distribution
	
	f = open(cosine_path,"w")
	# cos of angle from sat to Host		host/partner_massratio		eccentricity c/a(for binary) 	eccentricity b/a(for binary)	eccentricity c/b (for binary)	separation distance		dist of satellite to host		distance from satellite to host/dsep	sat mass
	f.write("cos_angle" +  "\t" + "host/partner_massratio" + "\t" + "c/a" + "\t" + "b/a" + "\t" + "c/b" + "\t" + "dsep" + "\t" + "dS" + "\t" + "dS/dsep" + "\t" + "mS" + "\n")
	for j in range(len(LGID)):
		if j % 2 == 0:
			pairno = int(j/2e0) #unpack the data stored in "_pop.dat" in pairs (pairs were written sequentially)
			print ("pair number", pairno)
			small_halo_path = os.path.join(os.getcwd(), folder_control, "LGfakeBIN" + str(pairno) + "_small_halos_feats.dat")
			SMALLx, SMALLy, SMALLz, SMALLm = loadtxt(small_halo_path, skiprows = 0, unpack=True)#get small halo positions
			
			SMALL = [SMALLx, SMALLy, SMALLz, SMALLm]
	
			mA, mB = LGm[j], LGm[j+1] #host mass
			bA, bB = LGb[j], LGb[j+1] #host eccentricity; b is b/a
			cA, cB = LGc[j], LGc[j+1] #host eccentricity; c is c/a

			Sepdis = sepDistlist[int(j/2)]
			hostApos = np.array([0, 0, 0])
			hostBpos = np.array([Sepdis, 0, 0])

			cosA, cosB = cos_satellite_angle_OL(mA, mB, hostApos, hostBpos, Sepdis, SMALL)

			# note: cosA[i] = [cosangle, dist to host, mass]
			for i in range(len(cosA)):
				# file line: cos_angle	hostmassratio	eccentricity b/c 	eccentricity b/a	satellite mass/host mass]
				f.write(str(cosA[i][0]) + "\t" + str(mA/mB) + "\t" + str(cA) + "\t" + str(bA) + "\t" + str(cA/bA) + "\t" + str(Sepdis) + "\t" + str(cosA[i][1]) + "\t" + str(cosA[i][1]/Sepdis) +"\t" + str(cosA[i][2]) + "\n")
			for i in range(len(cosB)):
				f.write(str(cosB[i][0]) + "\t" + str(mB/mA) + "\t" + str(cB) + "\t" + str(bB) + "\t" + str(cB/bB) + "\t" + str(Sepdis) + "\t" + str(cosB[i][1]) + "\t" + str(cosB[i][1]/Sepdis) +"\t" + str(cosB[i][2]) + "\n")
	f.close()

def cosdistn(LGID, LGm, LGx, LGy, LGz, LGRvir, LGb, LGc, cosine_path):
	#plotting histogram of the cosine angle distribution
	
	f = open(cosine_path,"w")
	# cos of angle from sat to Host		host/partner_massratio		eccentricity c/a(for binary) 	eccentricity b/a(for binary)	eccentricity c/b (for binary)	separation distance		dist of satellite to host		distance from satellite to host/dsep	sat mass
	f.write("cos_angle" +  "\t" + "host/partner_massratio" + "\t" + "c/a" + "\t" + "b/a" + "\t" + "c/b" + "\t" + "dsep" + "\t" + "mS"  + "\t" + "xS" + "\t" + "yS" + "\t" + "zS"+ "\t" + "dS/dsep" + "\n")
	
	for j in range(len(LGID)):
		if j % 2 == 0:
			pairno = int(j/2e0) #unpack the data stored in "_pop.dat" in pairs (pairs were written sequentially)
			print ("pair number", pairno)
			if pairno not in neglist:
				small_halo_path = os.path.join(os.getcwd(), folder_smallhalo, str(pairno) + "_small_halos_feats.dat")
				SMALLm, SMALLx, SMALLy, SMALLz = loadtxt(small_halo_path, skiprows = 1, usecols = (masscol, xcol, ycol, zcol), unpack=True)#get small halo positions
				SMALL = [SMALLm, SMALLx, SMALLy, SMALLz]
			
				mA, mB = LGm[j], LGm[j+1] #host mass
				bA, bB = LGb[j], LGb[j+1] #host eccentricity; b is b/a
				cA, cB = LGc[j], LGc[j+1] #host eccentricity; c is c/a

				hostApos = np.array([LGx[j], LGy[j], LGz[j]])
				hostBpos = np.array([LGx[j+1], LGy[j+1], LGz[j+1]])
				
				Sepdis = coorddistance(hostApos, hostBpos)
				cosA, cosB = cos_satellite_angle(mA, mB, hostApos, hostBpos, Sepdis, SMALL)
		
				# note: cosA[i] = [cosangle, mass of satellite, dist to host, satx, saty, satz]
				for i in range(len(cosA)):
				# file line: cos_angle	hostmassratio	eccentricity b/c 	eccentricity b/a	eccentricity c/b  sep DIs		mass sat	d sat/ dsep
					f.write(str(cosA[i][0]) + "\t" + str(mA/mB) + "\t" + str(cA) + "\t" + str(bA) + "\t" + str(cA/bA) + "\t" + str(Sepdis) + "\t" + str(cosA[i][1]) + "\t" + str(cosA[i][3]) + "\t" + str(cosA[i][4]) + "\t" + str(cosA[i][5]) + "\t" + str(cosA[i][2]/Sepdis) + "\n")
				for i in range(len(cosB)):
					f.write(str(cosB[i][0]) + "\t" + str(mB/mA) + "\t" + str(cB) + "\t" + str(bB) + "\t" + str(cB/bB) + "\t" + str(Sepdis) + "\t" + str(cosB[i][1]) + "\t" + str(cosB[i][3]) + "\t" + str(cosB[i][4]) + "\t" + str(cosB[i][5]) + "\t" + str(cosB[i][2]/Sepdis) + "\n")
	f.close()
		
if __name__ == "__main__":
	Rockstarfilelist0 = []
	Rockstar_files = os.listdir(Rockstar_hlist_address)
	
	for i in range(len(Rockstar_files)):
		filename = Rockstar_files[i]
		if filename.endswith("list"):
			Rockstarfilelist0.append(filename)
	Rockstarfilelist0.sort(key = natural_keys)
	filename_high_res = Rockstarfilelist0[::-1][0]
	
	# copy hlist file to home address
	snapnum = filename_high_res[6 : 12]
	datfoldr = "snapshot" + snapnum
		
	filepath_catalog = Rockstar_hlist_address + filename_high_res
	newfilepath = datfoldr + "/" + filename_high_res
#	copyfile(filepath_catalog, newfilepath)
	os.chdir("./" + datfoldr)

	neglist = loadtxt("negfile.dat")

	hrLGcathaloID, hrLGm, hrLGx, hrLGy, hrLGz, hrLGRvir, hrLGb, hrLGc = col_reader5(high_res_trueLG_pop, IDcol, masscol, xcol, ycol, zcol, Rvircol, bcol, ccol)
	
	cosdistn(hrLGcathaloID, hrLGm, hrLGx, hrLGy, hrLGz, hrLGRvir, hrLGb, hrLGc, hrcosine_distn_file) #for real signals
	
###############################################################
	hrfakeLGID, hrfakelgm, hrfakelgx, hrfakelgy, hrfakelgz, hrfakelgRvir, hrfakelgb, hrfakelgc = col_reader5(filename_hr_fakeBIN_pop, IDcol, masscol, xcol, ycol, zcol, Rvircol, bcol, ccol)
	sepDistlist = loadtxt('fakepairlist.dat', dtype = float, usecols = (2,), unpack = True)
##	sys.exit()
	cosdistn_fake(hrfakeLGID, hrfakelgm, hrfakelgx, hrfakelgy, hrfakelgz, hrfakelgRvir, hrfakelgb, hrfakelgc, sepDistlist, hrcosine_distn_file_fake, folder_control) #produces cosine distn file


