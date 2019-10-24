import os
from numpy import loadtxt
import numpy as np
from constant import *
import re
import sys
from numpy import sqrt

#GENERATE MARKOV RANDOM DISTRIBUTION

def markovN2(N):
	q_list = []
	for i in range(10000):
		s = np.random.uniform(-1e0, 1e0, N) #
		NA = np.where(np.logical_and(s < 1.0, s > 0))
		NB = np.where(np.logical_and(s > -1.0, s < 0))
		q_list.append(len(NA[0])/len(NB[0]))
		
	stdAB = np.std(q_list)
	meanAB = np.mean(q_list)
	return meanAB, stdAB

def markovN(N):
	q_list = []
	for i in range(10000):
		s = np.random.uniform(-1e0, 1e0, N) #
		NA = np.where(np.logical_and(s < 1.0, s > 0.8))
		NB = np.where(np.logical_and(s > -1.0, s < -0.8))
		q_list.append(len(NA[0])/len(NB[0]))
		
	stdAB = np.std(q_list)
	meanAB = np.mean(q_list)
	return meanAB, stdAB

def markovNabs(N):
	meanlist = np.array([])
	for i in range(10000):
		s = np.random.uniform(0, 1e0, N) #
		mean_s = np.mean(s)
		np.append(meanlist, mean_s)
	meanabs = np.mean(meanlist)
	stdabs = np.std(meanlist)
	return meanabs, stdabs
		
def markov(N):
	medians_lst = []
	for i in range(10000):
		s = np.random.uniform(-1e0, 1e0, N)
		median = np.median(s)
		medians_lst.append(median)
	return (np.median(medians_lst), np.std(medians_lst))
		
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

####################### DISTANCE CALCULATORS AND COORDINATE WRAPPER ###################################
#def middlecoordmpc2(coord1, coord2):
#	return np.mod(((np.mod(np.mod(coord1+coord2,dim), dim))/ 2e0),dim)


def lawofcos(coordS, coordA, coordB):
	# explanation to the function law of cosine:
	# cos A = (b^2 + c^2 - a^2)/2bc
	# A is the angle of the line from satellite to the host closest to it to the line connecting two hosts	
	
	b_sq = coorddistance_sqd(coordS, coordA)
	c_sq = coorddistance_sqd(coordB, coordA)
	a_sq = coorddistance_sqd(coordS, coordB) #B needs to be the host further to S (satellite) in this function
	return (b_sq + c_sq - a_sq) / (2e0 * sqrt(b_sq) * sqrt(c_sq) )

def lawofcos_arr(coordS, coordA, coordB):
	# lawofcos for array
	b_sq = coorddistance_sqd_arr(coordS, coordA)
	c_sq = coorddistance_sqd_arr(coordB, coordA)
	a_sq = coorddistance_sqd_arr(coordS, coordB) #B needs to be the host further to S (satellite) in this function
	return (b_sq + c_sq - a_sq) / (2e0 * sqrt(b_sq) * sqrt(c_sq) )
		
def middlecoordmpc(coord1, coord2):
	md = []
	for (a, b) in zip(coord1, coord2):
		delta = abs(b - a)
		if delta > dim - delta: #if say a is 50, b is 9, they need to be wrapped around, delta is then 41
			md.append( coordwrap2(min(a,b) - (dim - delta) / 2e0) ) #then the midpoint should be (i.e. b) - half the wrapped distance 50 64 0 9
		else:
			md.append((a + b) / 2e0)
	return md

def coorddistance_sqd(coord1, coord2):
	return np.sum((np.mod(np.fabs(coord1-coord2)-dim/2e0, dim)-dim/2e0) ** 2e0, axis = 0)

def coorddistance_sqd_arr(coord1, coord2):
	return np.sum((np.mod(np.fabs(coord1-coord2)-dim/2e0, dim)-dim/2e0) ** 2e0, axis = 1)
				
def coorddistance(coord1, coord2):
	return (np.sum((np.mod(np.fabs(coord1-coord2)-dim/2e0, dim)-dim/2e0) ** 2e0, axis = 0))**5e-1

def coorddistance_arr(coord1, coord2):
	return (np.sum((np.mod(np.fabs(coord1-coord2)-dim/2e0, dim)-dim/2e0) ** 2e0, axis = 1))**5e-1
	
def coordsubstract(coord1, coord2):
	diff = []
	for (a, b) in zip(coord1, coord2):
		delta = a - b
		if delta > 0 and delta > dim - delta: #a is 64, b is 1
			delta = - dim + delta
		if delta < 0 and abs(delta) > dim + delta: #a is 1, b is 64
			delta = dim + delta
		diff.append(delta)
	return diff

def coordsubstract_arr(coord1, coord2):
	return (np.mod(coord1-coord2-dim/2e0, dim)-dim/2e0)

def coordcompare(a, b):
	#return 1 if a>b, or a is to the right of b in a wrapped box
	delta = a - b
	if delta > 0:
		if delta < dim - delta: #a is 64, b is 1
			return 1
		else:
			return 0
	if delta < 0:
		if abs(delta) < dim + delta: #a is 1, b is 64
			return 0
		else:
			return 1

def coordingrid(coord, grid):
#	if 3d coord is in a box defined by 3 x 2 matrix grid
	counter = 0
	for l in range(3): #for x y z					
		if grid[l][0] < grid[l][1]: #not wrap condition for coordwrap2 function: boxmin < boxmax
			if coord[l] > grid[l][0] and coord[l] < grid[l][1]: #larger than boxmin and smaller than boxmax
				counter += 1
			else: 
				return 0
		else: #wrap condition applies: boxmax < boxmin
			if coord[l] > grid[l][0] or coord[l] < grid[l][1]: #larger than boxmin or smaller than boxmax
				counter += 1
			else:
				return 0				
	if counter == 3: # coord is in the box
		return 1

def coordwrap(coord1, coord2):
	#wrap coordinate coord2 (which is usually at the far end of box) to be nearer to coord1 (near 0Mpc)
	marker = [] #stays zero if they don't have to be wrapped
	for j in range(3):
		delta = abs(coord1[j] - coord2[j])
		if delta > dim - delta: #they satisfy the wrap condition
			if coord1[j] > coord2[j]:
				marker.append(1) #coord1 is around 64mpc, after wrap it's a small negative
				coord1[j] = coord1[j] - dim
			else:
				marker.append(2) #coord2 is around 64mpc
				coord2[j] = coord2[j] - dim
		else:
			marker.append(0) #in j direction no wrap is required
	return coord1, coord2, marker
	
def coordwrap2(cd):
	#wrap out of simulation box coord to be inside
	return np.mod(cd, dim)

def coordwrap3(hostApos, hostBpos, marker, pos):
	#wrap a point to be near 0
	for j in range(3):
		if marker[j] == 1: #A is wrapped to be near B
			delta = abs(pos[j] - hostBpos[j])
			if delta > dim - delta:
				if pos[j] > hostBpos[j]:
					pos[j] = pos[j] - dim
		if marker[j] == 2: #B is wrapped to be near A
			delta = abs(pos[j] - hostApos[j])
			if delta > dim - delta:
				if pos[j] > hostApos[j]:
					pos[j] = pos[j] - dim
		if marker[j] == 0: #A and B are not wrapped
			delta = abs(pos[j] - hostApos[j])
			if delta > dim - delta:
				if pos[j] > hostApos[j]:
					pos[j] = pos[j] - dim
				elif pos[j] < hostApos[j]:		
					pos[j] = pos[j] + dim
	return pos

def shiftpos(coord1, coord2):
#	shifting coord1 by coord2
	return coord1 - coord2 + 32e0

def shiftpos2(coord1, coord2, dist):
#	shifting coord1 by coord2 then by a distance in one direction
	return coord1 - coord2 + 32e0 + np.array([dist,0,0])

def threeDvel(vx, vy, vz):
	return (vx ** 2e0 + vy ** 2e0 + vz ** 2e0) ** .5e0
################################# FILE LOADERS ###################################

def xyz_converter(x0ar, y0ar, z0ar):
	return x0ar/1e3, y0ar/1e3, z0ar/1e3

def xy_converter(x0ar, y0ar):
	return x0ar/1e3, y0ar/1e3
		
#def xyz_reader(filename, xcol, ycol, zcol):
#	xar, yar, zar = loadtxt(filename, skiprows = 1, usecols = (xcol, ycol, zcol), unpack=True)
#	return xar/1e3, yar/1e3, zar/1e3

def xyz_reader2(filename, xcol, ycol, zcol):
	xar, yar, zar = loadtxt(filename, skiprows = 1, usecols = (xcol, ycol, zcol), unpack=True)
	return xar, yar, zar

def col_reader0(filename, masscol, xcol, ycol, zcol):
	m, x0ar, y0ar, z0ar = loadtxt(filename, skiprows = 1, usecols = (masscol, xcol, ycol, zcol), unpack=True)
	return m, x0ar, y0ar, z0ar

def col_reader(filename, IDcol, masscol, xcol, ycol, zcol, Rvircol):
	IDar = loadtxt(filename, skiprows = 1, usecols = [IDcol], unpack=True, dtype='i8') 
	m, x0ar, y0ar, z0ar, rvir = loadtxt(filename, skiprows = 1, usecols = (masscol, xcol, ycol, zcol, Rvircol), unpack=True)
	return IDar, m, x0ar, y0ar, z0ar, rvir 

def col_reader2(filename, IDcol, masscol, xcol, ycol, zcol):
	IDar = loadtxt(filename, skiprows = 1, usecols = [IDcol], unpack=True, dtype='i8') 
	m, x0ar, y0ar, z0ar = loadtxt(filename, skiprows = 1, usecols = (masscol, xcol, ycol, zcol), unpack=True)
	return IDar, m, x0ar, y0ar, z0ar

def col_reader3(filename, IDcol, masscol, xcol, ycol, zcol, Rvircol):
	IDar = loadtxt(filename, skiprows = 1, usecols = [IDcol], unpack=True, dtype='i8') 
	m, x0ar, y0ar, z0ar, rvir = loadtxt(filename, skiprows = 1, usecols = (masscol, xcol, ycol, zcol, Rvircol), unpack=True)    
	return IDar, m, x0ar/1e3, y0ar/1e3, z0ar/1e3, rvir
		
def col_reader4(filename, IDcol, hostIDcol, masscol, xcol, ycol, zcol):
	IDar = loadtxt(filename, skiprows = 1, usecols = [IDcol], unpack=True, dtype='i8') #dtype='i8' makes reading long string like the ID number, possible	
	HOSTIDar = loadtxt(filename, skiprows = 1, usecols = [hostIDcol], dtype='str') 
	HOSTIDlist = []
	for hostID in HOSTIDar:
		if hostID == '0':
			HOSTIDlist.append(0)
		else: 
			HOSTIDlist.append(long(hostID, 10))
			
	m, x0ar, y0ar, z0ar = loadtxt(filename, skiprows = 1, usecols = (masscol, xcol, ycol, zcol), unpack=True) 
	return IDar, HOSTIDlist, m, x0ar/1e3, y0ar/1e3, z0ar/1e3

def col_reader5(filename, IDcol, masscol, xcol, ycol, zcol, Rvircol, bcol, ccol):
	IDar = loadtxt(filename, skiprows = 1, usecols = [IDcol], unpack=True, dtype='i8') #dtype='i8' makes reading long string like the ID number, possible
	m, x0ar, y0ar, z0ar, rvir, b, c = loadtxt(filename, skiprows = 1, usecols = (masscol, xcol, ycol, zcol, Rvircol, bcol, ccol), unpack=True) 
	
	return IDar, m, x0ar, y0ar, z0ar, rvir, b, c

def col_reader6(filename, IDcol, descIDcol, masscol, xcol, ycol, zcol, vxcol, vycol, vzcol, Rvircol, bcol, ccol):
	ID, descID = loadtxt(filename, skiprows = 1, usecols = [IDcol, descIDcol], unpack=True, dtype='i8') 
	m, x, y, z, vx, vy, vz, rvir, b, c = loadtxt(filename, skiprows = 1, usecols = (masscol, xcol, ycol, zcol, vxcol, vycol, vzcol, Rvircol, bcol, ccol), unpack=True) 
	return ID, descID, m, x, y, z, vx, vy, vz, rvir, b, c

def col_reader7(filename, IDcol, descIDcol, masscol, xcol, ycol, zcol, vxcol, vycol, vzcol):
	ID, descID = loadtxt(filename, skiprows = 1, usecols = [IDcol, descIDcol], unpack=True, dtype='i8') 
	m, x, y, z, vx, vy, vz = loadtxt(filename, skiprows = 1, usecols = (masscol, xcol, ycol, zcol, vxcol, vycol, vzcol), unpack=True) 
	return ID, descID, m, x, y, z, vx, vy, vz

def col_reader8(filename, IDcol, descIDcol, masscol, xcol, ycol, zcol):
	ID, descID = loadtxt(filename, skiprows = 1, usecols = [IDcol, descIDcol], unpack=True, dtype='i8') 
	m, x, y, z = loadtxt(filename, skiprows = 1, usecols = (masscol, xcol, ycol, zcol), unpack=True) 
	return ID, descID, m, x, y, z
			
def readindex(filename):
	#read indices of the binary LG from their index contained in the datafile lowresLG.dat and highresLG.dat
	jlst, klst = loadtxt(filename, skiprows=1, usecols = (0,1), unpack=True, dtype = int)
	dlst = loadtxt(filename, skiprows=1, usecols = (2,), unpack=True)
	return jlst, klst, dlst

def readindex2(filename):
	#read indices of the binary LG from their index contained in the datafile lowresLG.dat and highresLG.dat
	jlst = loadtxt(filename, skiprows=1, unpack=True, dtype = int)
	return jlst
		
def readpopinfo(filteredhalofile, trueLGhalofile, jlst, klst):
	# Read data from filtered halo file
	f = open(filteredhalofile,"r")
	lines = f.readlines()
	f.close()
	f2 = open(trueLGhalofile,"w")
	f2.write(lines[0])# write header
#	for i in range(65):
#		print(lines[i])
	for h in range(len(jlst)):
		f2.write(lines[jlst[h]+1]) #+headerlength because of header lines
		f2.write(lines[klst[h]+1])
	f2.close()

def readpopinfo2(filteredhalofile, trueLGhalofile, jlst):
	print(len(jlst))
	# Read data from filtered .AHF_halos file
	f = open(filteredhalofile,"r")
	lines = f.readlines()
	f.close()
	ID = loadtxt(filteredhalofile, skiprows = 1, usecols = (IDcol,), unpack=True, dtype='i8')	
	f2 = open(trueLGhalofile,"w")
	f2.write(lines[0])# write header
	for h in range(len(jlst)):
		f2.write(lines[int(jlst[h])+1]) #+1 because of header line
	f2.close()
	print("AHF entries of the " + str(len(jlst)) + " true isolated haloes are saved under: " + trueLGhalofile)

def massratio(m):
	massratiolst = []
	for j in range(len(m)):
		if j % 2 == 0:
			massratio = m[j]/m[j+1]
			if massratio > 1:
				massratio = m[j+1]/m[j]
			massratiolst.append(massratio)
	return massratiolst

def masssum(m):
	masssumlst = []
	for j in range(len(m)):
		if j % 2 == 0:
			masssum = m[j] + m[j+1]
			masssumlst.append(masssum)
	return masssumlst
	
def distanceLG(coord_x, coord_y, coord_z):
	distancelst = []
	for j in range(len(coord_x)):
		if j % 2 == 0:			
			distance = coorddistance([coord_x[j+1],coord_y[j+1],coord_z[j+1]],[coord_x[j],coord_y[j],coord_z[j]])
			distancelst.append(distance)
	return distancels
	
def paperaddress(pngname):
	return os.path.join(os.getcwd(), "paper", pngname)
	
def gapfilteredfilewriter(LGBx):
	f2 = open(high_res_trueLG_pop,"r") #read header
	lines = f2.readlines()
	f2.close()
	
	f = open(hrtrueLG_wosubhalo_biggap_filtered,"w")
	f.write(lines[0])
	
	neglist = loadtxt("negfile.dat")
	
	for i in range(len(LGBx)):
		if i % 2==0:
			pairno = int(i/2e0)
			if pairno not in neglist:
				f.write(lines[i+1])
				f.write(lines[i+2])
	f.close()
	
	
def propagate(A, B, sigma_A, sigma_B):
	#propagate errors of A/B
	
	return( A/B * sqrt((sigma_A/A)**2e0 + (sigma_B/B)**2e0) )

def propagate2(A, B, sigma_A, sigma_B):
	#propagate errors of A/B
	R=A/B
	P=(A+sigma_A) / (B-sigma_B)
	Q=(A-sigma_A) / (B+sigma_B)
	p= abs(R-P)
	q=abs(R-Q)
	r=max([p,q])

	return (r)
	
