'''
Produce a list of haloes that are isolated in space from other LG mass haloes and OM mass haloes, used for the purpose of measuring overlapping effect. 
'''
import sys
import os
from shutil import copyfile

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import sqrt, random
from constant import *
from helpers import *
from checknosubhalo import *
from subhaloFinderIndexed import *
from overlappingplot import *
from histplot import *
from findprince import findprince

def random3DRotMat():
	[x,y,z] = np.random.rand(3,1000)
	A = findprince(x,y,z)
	for v in A:
		vlength = (v[0] ** 2e0 + v[1] ** 2e0 + v[2] ** 2e0) ** .5e0
		for i in range(3):
			v[i] = v[i] /vlength
	return A
		
def rotate(x, y, z, v):	
	# v is vector to be rotated by, is a normalized vector
	
	N = len(x)
	nx = np.zeros(N)
	ny = np.zeros(N)
	nz = np.zeros(N)
	
	for i in range(N):
#		[nx[i], ny[i], nz[i]] = dot(np.array([x[i],y[i],z[i]]), v).tolist()
		nx[i] = v[0][0]*x[i] + v[1][0]*y[i] + v[2][0]*z[i]
		ny[i] = v[0][1]*x[i] + v[1][1]*y[i] + v[2][1]*z[i]
		nz[i] = v[0][2]*x[i] + v[1][2]*y[i] + v[2][2]*z[i]
	
	return np.array(nx), np.array(ny), np.array(nz)
	
def SimDist(Mi, MH, bi, bH, ci, cH): #similarity distance between two halos
	d = ((MH - Mi)/max(MH, Mi)) ** 2e0 + ((bH - bi)/max(bH, bi)) ** 2e0 + ((cH - ci)/max(cH, ci)) ** 2e0
#	print ((MH - Mi)/max(MH, Mi)) ** 2e0 #only selecting mass -> not a good effect
	return d

def truISOLG_finder():

	# m, x, y, z are the isolated candidates, omx, omy, omz are the overmassive haloes
	# we want to make sure the haloes among the candidates which stand alone with respect both to other candidates 
	# as well as to the overmassive ones are chosen, fixed distance is used to get a feeling for the sample
	filename = isolated_lg_file_at_a_Dsep
	
	f3 = open(filename, "w")
	f3.write("index of isolated halo/ line in file _isolated_lg_halo_nosub.dat" +  "\t" + "index of the binary halo in the sample associated with this isolation radius" + "\n")
	f3.close()
#	
	sepDis_LGBj_list = []
	
	for j in range(len(Tlgbx)):
		if j % 2 == 0:
			coordlgbA = np.array([Tlgbx[j], Tlgby[j], Tlgbz[j]])
			coordlgbB = np.array([Tlgbx[j + 1], Tlgby[j + 1], Tlgbz[j + 1]])	
			massA = Tlgbm[j]
			massB = Tlgbm[j + 1]
			sepDis = coorddistance(coordlgbA, coordlgbB)
			sepDis_LGBj_list.append((j / 2, sepDis, [massA, massB]))
			print(sepDis)
			
	sepDis_LGBj_list.sort(key = lambda tup: tup[1], reverse=True) # sort by sep distance from large to small (isolation criteria hard to easy)
	
	ISO_alreadyused = [] #stores true ISO index which has been used already
	
	for j in range(len(sepDis_LGBj_list)):
		LGB_index = sepDis_LGBj_list[j][0]
		sepDis = sepDis_LGBj_list[j][1]
		massA = sepDis_LGBj_list[j][2][0]
		massB = sepDis_LGBj_list[j][2][1]		
		isosearcadius = sepDis * 1.5e0
		
		pos = 0 
		for i in range(len(isoLGx)):
			if isoLGRvir[i] / 1e3 < sepDis / 2e0:
				if i not in ISO_alreadyused:
					coord = np.array([isoLGx[i], isoLGy[i], isoLGz[i]])
					isomass = isoLGm[i]

	#				simd = 1e30 #just need to be some large number to start the loop
	#				simd_new = SimDist(isomass, lgbm[i], b[j], lgbb[i], c[j], lgbc[i])
	#				if simd_new < simd: #similarity distance in terms of mass, b and c
	#				simd = simd_new
	#				real = i
				
					boxmin = coordwrap2(coord - isosearcadius) #within a box extended by 1.5Mpc * 1.5 in tee direction 
					boxmax = coordwrap2(coord + isosearcadius) #with the jth halo at its center
					grid = np.transpose(np.array([boxmin, boxmax]))

					omm = np.ones(len(omx)) * 1e13 # as long as it's larger than the LGM_max is fine
					OTHERcoordx = np.concatenate((lgx, omx))
					OTHERcoordy = np.concatenate((lgy, omy))
					OTHERcoordz = np.concatenate((lgz, omz))
					OTHERm = np.concatenate((lgm, omm))
			
					OTHER = np.transpose(np.array([OTHERcoordx, OTHERcoordy, OTHERcoordz]))
				
					index_inbox = np.where( np.ndarray.all(np.mod(OTHER - boxmin - dim / 2e0, dim) - dim / 2e0 > 0, axis = 1) & np.ndarray.all(np.mod(OTHER - boxmax - dim / 2e0, dim) - dim / 2e0 < 0, axis = 1) )
					index = np.where( np.logical_and(coorddistance_sqd_arr(OTHER[index_inbox[0]], np.array([coord])) < isosearcadius ** 2e0,  coorddistance_sqd_arr(OTHER[index_inbox[0]], np.array([coord])) != 0) )
	#				if j == 1 and i == 1949:
	#					print(len(OTHER), len(omx), len(lgx))
	#					print(sepDis, OTHERm[index_inbox[0][index[0]]], isomass, OTHER[index_inbox[0][index[0]]], coord)
					
					neg = 0
					for k in index_inbox[0][index[0]]: #of all second halos (of LG mass over overmassive) insdie the sphere of isosearcadius
						if OTHERm[k] >= min(massA, massB) or OTHERm[k] >= isomass:
							neg = 1
							break # halo i is not isolated at distance demanded by LGB pair LGB_index

					if neg == 0: # halo i is isolated at LGB_index's dsep
						f3 = open(filename, "a")
	#					if i == 1949:
	#						print('something')
	#						for k in index_inbox[0][index[0]]:
	#							print(coorddistance_sqd_arr(OTHER[k], np.array([coord])), OTHERm[k], m[i])
	#						sys.exit()
						if abs(isomass - massA) > abs(isomass - massB):
							f3.write(str(i) + "\t" + str(int(LGB_index)) + "\t" + str(isomass) + "\t" + str(sepDis) + "\t" + str(massB) + "\n")
						else:
							f3.write(str(i) + "\t" + str(int(LGB_index)) + "\t" + str(isomass) + "\t" + str(sepDis) + "\t" + str(massA) + "\n")
						
						ISO_alreadyused.append(i)
						f3.close()
						pos = 1 # there exist isolated halo candidates for isolation level at LGB_index's dsep

		if pos == 0:
			f3 = open(filename, "a")
			f3.write(str(-1) + "\t" + str(int(LGB_index)) + "\t" + str(0) + "\t" + str(sepDis) + "\t" + str(0) + "\n")
			f3.close()
	
def fakebin_finder():
	sepDis_list_orderedbyIndex = []
	for j in range(len(Tlgbx)):
		if j % 2 == 0:
			coordlgbA = np.array([Tlgbx[j], Tlgby[j], Tlgbz[j]])
			coordlgbB = np.array([Tlgbx[j + 1], Tlgby[j + 1], Tlgbz[j + 1]])	
			
			sepDis = coorddistance(coordlgbA, coordlgbB)
			sepDis_list_orderedbyIndex.append(sepDis)
			print(sepDis)
	
	iso_candidate, sample_halo_index = loadtxt(isolated_lg_file_at_a_Dsep, dtype = int, usecols = (0,1), unpack = True, skiprows = 1)
	iso_mass, dsep_list, host_mass = loadtxt(isolated_lg_file_at_a_Dsep, dtype = float, usecols = (2,3,4), unpack = True, skiprows = 1)
	
	fake_pair_list = []
	isoHalo_Alreadytaken = [-1]
	
	for j in range(len(Tlgbx)):
		if j % 2 == 0:
			pairlocat = np.where(sample_halo_index == j//2)
			dsep = dsep_list[pairlocat[0][0]]
			
			rotate_tag = '0' #no rotation
			
			dist = 1e30
			for i in range(pairlocat[0][0]): # all the isos which can be candidate in fake pair for real pair j/2 
				if iso_candidate[i] not in isoHalo_Alreadytaken:
					dist_new = (Tlgbm[j] - iso_mass[int(i)]) ** 2e0
					if dist_new < dist:
						closestinmasstoA = i
						dist = dist_new
			try:
				fakeA = closestinmasstoA
				
			except NameError:
				dist = 1e30
				for i in range(pairlocat[0][0]):
					dist_new = (Tlgbm[j] - iso_mass[int(i)]) ** 2e0
					if dist_new < dist:
						closestinmasstoA = i
						dist = dist_new			
				fakeA = closestinmasstoA
				rotate_tag = 'A'

			dist = 1e30
			for i in range(pairlocat[0][0]):
				if iso_candidate[i] not in isoHalo_Alreadytaken:
					dist_new = (Tlgbm[j+1] - iso_mass[int(i)]) ** 2e0
					if dist_new < dist:
						closestinmasstoB = i
						dist = dist_new
			try:
				fakeB = closestinmasstoB
				
			except NameError:
				dist = 1e30
				for i in range(pairlocat[0][0]):
					dist_new = (Tlgbm[j+1] - iso_mass[int(i)]) ** 2e0
					if dist_new < dist:
						closestinmasstoB = i
						dist = dist_new			
				fakeB = closestinmasstoB
				if rotate_tag == '0':
					rotate_tag = 'B'
				else:
					rotate_tag = 'AB'

			if rotate_tag == '0':
				fake_pair_list.append([iso_candidate[fakeA], iso_candidate[fakeB], dsep, 0, Tlgbm[j], Tlgbm[j+1], iso_mass[fakeA], iso_mass[fakeB]])
			elif rotate_tag == 'A':
				fake_pair_list.append([iso_candidate[fakeA], iso_candidate[fakeB], dsep, 1, Tlgbm[j], Tlgbm[j+1], iso_mass[fakeA], iso_mass[fakeB]])			
			elif rotate_tag == 'B':
				fake_pair_list.append([iso_candidate[fakeA], iso_candidate[fakeB], dsep, 2, Tlgbm[j], Tlgbm[j+1], iso_mass[fakeA], iso_mass[fakeB]])
			elif rotate_tag == 'AB':
				fake_pair_list.append([iso_candidate[fakeA], iso_candidate[fakeB], dsep, 3, Tlgbm[j], Tlgbm[j+1], iso_mass[fakeA], iso_mass[fakeB]])
										
#			if pairno < 10:
#				print(fake_pair_list)
	fake_pair_list = np.transpose(np.array(fake_pair_list))
	DAT = np.column_stack((fake_pair_list[0], fake_pair_list[1], fake_pair_list[2], fake_pair_list[3], fake_pair_list[4], fake_pair_list[5], fake_pair_list[6], fake_pair_list[7]))
	np.savetxt('fakepairlist.dat', DAT, delimiter=" ")

def truISOLG_sat_finder():
	folder = "smallhalos_highres_iso"
	try:
		os.makedirs(folder)
	except OSError:
		pass
	
	fakeA_col, fakeB_col, tag_col = loadtxt('fakepairlist.dat', usecols = (0,1,3), unpack = True)
	dsep_list = loadtxt('fakepairlist.dat', dtype = float, usecols = (2,), unpack = True)
		
	for j in range(len(fakeA_col)):
		print(j)
		fakeA_ind = int(fakeA_col[j])
		fakeB_ind = int(fakeB_col[j])
		rotation_tag = int(tag_col[j])
		dsep = dsep_list[j]
		isosearcadius = dsep * 1.5e0
		print(dsep)
		coord_fakeA = [isoLGx[fakeA_ind], isoLGy[fakeA_ind], isoLGz[fakeA_ind]]
		coord_fakeB = [isoLGx[fakeB_ind], isoLGy[fakeB_ind], isoLGz[fakeB_ind]]
		
		s_fAB = []
		smass_fAB = []
		ind = [fakeA_ind, fakeB_ind]
		
		for i, coord in enumerate([coord_fakeA, coord_fakeB]):
			boxmin = coordwrap2(coord - isosearcadius) #within a box extended by 1.5Mpc * 1.5 in tee direction 
			boxmax = coordwrap2(coord + isosearcadius) #with the jth halo at its center
			grid = np.transpose(np.array([boxmin, boxmax]))

			SMALL = np.transpose(np.array([totx, toty, totz]))
		
			index_sat_inbox0 = np.where( np.ndarray.all(np.mod(SMALL - boxmin - dim / 2e0, dim) - dim/2e0 > 0, axis = 1) & np.ndarray.all(np.mod(SMALL - boxmax - dim/2e0, dim) - dim/2e0 < 0, axis = 1) )
			index_sat_inbox = index_sat_inbox0[0]
			index_sat0 = np.where( (coorddistance_sqd_arr(SMALL[index_sat_inbox], np.array([coord])) < isosearcadius**2e0) & (coorddistance_sqd_arr(SMALL[index_sat_inbox], coord) != 0 ) )
			index_sat = index_sat0[0]
			tot_inbox = index_sat_inbox[index_sat]
			smallcoord_f = SMALL[tot_inbox]
			smallmass_f = totm[tot_inbox]
			###### check point for isolation of the halo
			
			index_FORBIDDEN0 = np.where(totm[tot_inbox] > isoLGm[ind[i]])
			index_FORBIDDEN = index_FORBIDDEN0[0]
			if index_FORBIDDEN != []:
				print(dsep, ind[i],  totm[tot_inbox[index_FORBIDDEN]], isoLGm[ind[i]], SMALL[tot_inbox[index_FORBIDDEN]], coord)
				print("something seriously wrong! Halo not isolated")
				sys.exit()
				
			satpos0 = coordsubstract_arr(smallcoord_f, coord)
			satpos0 = np.transpose(satpos0)
			sx = satpos0[0]
			sy = satpos0[1]
			sz = satpos0[2]
			
			A_random = random3DRotMat()

			if rotation_tag == 1 and i == 0: #A rotate
				sx, sy, sz = rotate(sx, sy, sz, A_random)
							
			if rotation_tag == 2 and i == 1: #B rotate
				sx, sy, sz = rotate(sx, sy, sz, A_random)

			if rotation_tag == 3: #both rotate
				sx, sy, sz = rotate(sx, sy, sz, A_random)
			
			#shift B's satellite	
			if i == 1:
				sx += dsep
		
			s_fAB.append([sx, sy, sz])
			smass_fAB.append(smallmass_f)
			
		sx_fAB = np.concatenate((s_fAB[0][0], s_fAB[1][0]))
		sy_fAB = np.concatenate((s_fAB[0][1], s_fAB[1][1]))		
		sz_fAB = np.concatenate((s_fAB[0][2], s_fAB[1][2]))
		sm_fAB = np.concatenate((smass_fAB[0], smass_fAB[1]))
		
		small_halo_filename = "LGfakeBIN" + str(j) + "_small_halos_feats.dat"
		small_halo_path_fakeBIN = os.path.join(os.getcwd(), folder, small_halo_filename)
		
		DAT = np.column_stack((sx_fAB, sy_fAB, sz_fAB, sm_fAB))
		np.savetxt(small_halo_path_fakeBIN, DAT, delimiter=" ")

def fakepair_pop():
	f = open(isolated_lg_file_nosub,"r") #read header
	lines = f.readlines()
	lines = lines[1:]
	f.close()
	
	fakeA_col, fakeB_col = loadtxt('fakepairlist.dat', usecols = (0,1), unpack = True)
	f = open(filename_hr_fakeBIN_pop, "w")
	f.write(lines[0])
	print(len(fakeA_col))
	for i in range(len(fakeA_col)):
		f.write(lines[int(fakeA_col[i])])
		f.write(lines[int(fakeB_col[i])])		
	f.close()

def cumulativeFakeHostCumulative(cumtag):
	
	fakeA_col, fakeB_col, tag_col = loadtxt('fakepairlist.dat', usecols = (0,1,3), unpack = True)
	
	MRlist = []	
	for j in range(len(fakeA_col)):
		fakeA_ind = int(fakeA_col[j])
		fakeB_ind = int(fakeB_col[j])
		fakeA_m = isoLGm[fakeA_ind]
		fakeB_m = isoLGm[fakeB_ind]		
		MR = fakeA_m/fakeB_m
		if MR >1:
			MR = 1e0/MR
		MRlist.append(MR)
	MRmax = 1e0
	nbins = 40
	Massratiobins = np.arange(0, MRmax + MRmax/nbins, MRmax/nbins)
	
	histx, bin_edgesx = np.histogram(MRlist, bins=Massratiobins)
#	cumulative = np.cumsum(histx)

	cumulative = []
	for i in range(len(histx)):
		print(histx[i]/sum(histx))
		print('0', sum(histx[:i+1])/sum(histx))
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
	if cumtag ==1:
		ax.set_ylabel('$N_{pair}(<M_2/M_1)$',fontsize=20)
	else:
		ax.set_ylabel('$N_{pair}$',fontsize=20)
	plt.grid(True)
	if cumtag ==1:
		fig.savefig('Nfakehost_vs_Massratio_cumulative.jpg', dpi=500, bbox_inches='tight')
	else:
		fig.savefig('Nfakehost_vs_Massratio.jpg', dpi=500, bbox_inches='tight')
		
def ControlRepeatedOccurence_histplot():
	fig = plt.figure(figsize=(6,6))
	ax = fig.add_subplot(111)

	ax.set_xlabel("occurence",fontsize=14)
	ax.set_ylabel("N",fontsize=14)
	
	iso_ID, occur = loadtxt(filename_ControlRepeatedOccurence,unpack=True)
	plt.hist(occur, max(occur), facecolor='green')

	fig.savefig("ControlRepeatedOccurence.pdf", dpi=500)
	plt.close()

def ControlRepeatedOccurence_pair_histplot():
	fig = plt.figure(figsize=(6,6))
	ax = fig.add_subplot(111)

	ax.set_xlabel("occurence",fontsize=14)
	ax.set_ylabel("N",fontsize=14)
	
	occur = loadtxt(filename_ControlRepeatedOccurence_pair, usecols = (2,), unpack=True)
	plt.hist(occur, max(occur), facecolor='green')

	fig.savefig("ControlRepeatedOccurence_pair.pdf", dpi=500)
	plt.close()


if __name__ == "__main__":
	isoLGID, isoLGm, isoLGx, isoLGy, isoLGz, isoLGRvir, isob, isoc = col_reader5(isolated_lg_file_nosub, IDcol, masscol, xcol, ycol, zcol, Rvircol, bcol, ccol)
	print(len(isoLGID))
	cumulativeFakeHostCumulative(0)
	cumulativeFakeHostCumulative(1)	
	sys.exit()
##############################################################################################################################
############################################ begin of sample halo repeat checking ##############################################
##############################################################################################################################

############################side track: how repeated is the sample? Are there trios where A and B is a pair, as well as B and C
########################### Answer: yes. the pair from the halo search does not have duplicate, but individual haloes have duplicates 
########################### i.e. there are situations where AB, BC both qualify as LGB  #############################################

#	pairlist = []
#	duppairlist = []
#	for i in range(len(_jlst)):
#		pairlist.append((_jlst[i], _klst[i]))
#		pairlist.append((_klst[i], _jlst[i]))
#	print len(pairlist)/2
#	
#	for i in pairlist:
#		if pairlist.count(i)>1:
#			duppairlist.append(i)
#	print len(duppairlist)
# %no duplicate pairs

#	make a copy
#	jlst = []
#	klst = []
#	for i in _jlst:
#		jlst.append(i)
#	for i in _klst:
#		klst.append(i)
#				
#	duplist = []	
#	comblst = jlst + klst
#	for i in comblst:
#		if comblst.count(i)>1:
##			if i not in duplist:
#			duplist.append(i)
#	print len(duplist)

	#comblst = [1, 2, 3, 1, 3, 4, 2, 2]
	#duplist = [1, 2, 3, 1, 3, 2, 2]
	#jlst = [0, 2, 3, 1]
	#klst = [3, 4, 2, 2]
	
#	jilist = []
#	kilist = []
#	for i in duplist:
#		try:
#			ji = jlst.index(i)
#		except ValueError:
#			ji = -1
#		try:
#			ki = klst.index(i)
#		except ValueError:
#			ki = -1
#		
#		if ji >= 0:
#			jilist.append((i, ji))
#			jlst[ji] = 0
#		if ki >= 0:
#			kilist.append((i, ki))
#			klst[ki] = 0

#	filename_filtered = "high_res_massfiltered_lg.dat"
#	m, x0, y0, z0 = loadtxt(filename_filtered, skiprows = 1, usecols = (3, 5, 6, 7), unpack=True)	

#	tri_counter = 0
#	for (j, ji) in jilist: #jilist is a list of indices of things that have duplicate in either jlst or klst
#		for (k, ki) in kilist:
#			if j == k: #_jlst[ji] is the same as _klst[ki]
#				print "begin"
##				print j, ji
#				tri_counter +=1
#				print _jlst[ji], _klst[ji], sqrt(distlist[ji]), m[_jlst[ji]]/1e10, m[_klst[ji]]/1e10
#				print _jlst[ki], _klst[ki], sqrt(distlist[ki]), m[_jlst[ki]]/1e10, m[_klst[ki]]/1e10
#				coordA = [x0[_jlst[ji]]/1e3, y0[_jlst[ji]]/1e3, z0[_jlst[ji]]/1e3]
#				coordB = [x0[_klst[ji]]/1e3, y0[_klst[ji]]/1e3, z0[_klst[ji]]/1e3]
#				coordC = [x0[_jlst[ki]]/1e3, y0[_jlst[ki]]/1e3, z0[_jlst[ki]]/1e3]
#				c1 = middlecoordmpc(coordA, coordB)
#				c2 = middlecoordmpc(coordA, coordC)
#				print _jlst[ki], coorddistance(coordC, c2)
#				print _klst[ji], coorddistance(coordB, c1)
#	
#	print tri_counter

##############################################################################################################################
############################################ end of sample halo repeat checking ##############################################
##############################################################################################################################

#################################################### FIND ISOLATED HALOES ########################################################

	Rockstarfilelist0 = []
	Rockstar_files = os.listdir(Rockstar_hlist_address)
	
	for i in range(len(Rockstar_files)):
		filename0 = Rockstar_files[i]
		if filename0.endswith("list"):
			Rockstarfilelist0.append(filename0)
	Rockstarfilelist0.sort(key = natural_keys)
	filename_high_res = Rockstarfilelist0[::-1][0]
		
	# copy hlist file to home address
	snapnum = filename_high_res[6:12]
	datfoldr = "snapshot" + snapnum
		
	filepath_catalog = Rockstar_hlist_address + filename_high_res
	newfilepath = datfoldr + "/" + filename_high_res
#	copyfile(filepath_catalog, newfilepath)
	os.chdir("./" + datfoldr)
		
#	os.remove(filename_high_res)
	
	#step 1
	#all sample file of LG mass range haloes, call it file A
#	fA = open(filename_high_res_filtered_lg_nosub, "r")
#	lines0 = fA.readlines()
#	fA.close()
#	lines = lines0[1:]

#	fA2 = open(hrtrueLG_wosubhalo_biggap_filtered, "r") # with 4500 some LGB (some with large satellites)
#	lines2 = fA2.readlines()
#	fA2.close()
#	lines2 = lines2[1:]
#	
##	step 2
##	complimentary set of (A-B), do not contain subs, but contain isolated and non-isolated halos, called file C, C = A-B
#	fC = open(isolated_lg_file_nosub, "w")
#	fC.write(lines0[0])
#	for i in range(len(lines)):
#		if lines[i] not in lines2:
#			fC.write(lines[i])
#	fC.close()
	
#   step 3: check against each other and isolated halo to satisfy isolation condition 
	#unpack the no_sub file

	isoLGID, isoLGm, isoLGx, isoLGy, isoLGz, isoLGRvir, isob, isoc = col_reader5(isolated_lg_file_nosub, IDcol, masscol, xcol, ycol, zcol, Rvircol, bcol, ccol)
	print(len(isoLGID))
	
	#true binary info:
	TLGID, Tlgbm, Tlgbx, Tlgby, Tlgbz, TlgRvir, Tlgb, Tlgc = col_reader5(hrtrueLG_wosubhalo_biggap_filtered, IDcol, masscol, xcol, ycol, zcol, Rvircol, bcol, ccol)
	print(len(TLGID))
	
#	f = open('BinaryHostHalo_ID_mass_dsep.dat', "w")
#	f.write("halo ID" + "\t" + "mass" + "\t" + "d_sep" + "\n")
#	for i in range(len(TLGID)):
#		if i % 2 ==0:
#			coordlgbA = np.array([Tlgbx[i], Tlgby[i], Tlgbz[i]])
#			coordlgbB = np.array([Tlgbx[i+1], Tlgby[i+1], Tlgbz[i+1]])
#			sepDis = coorddistance(coordlgbA, coordlgbB)
#			f.write(str(TLGID[i]) + "\t" + str(Tlgbm[i]) + "\t" + str(sepDis) + "\n")
#			f.write(str(TLGID[i+1]) + "\t" + str(Tlgbm[i+1]) + "\t" + str(sepDis) + "\n")
#	f.close()
#	sys.exit()
	
	LGID, lgm, lgx, lgy, lgz, lgRvir, lgb, lgc = col_reader5(filename_high_res_filtered_lg_wsub, IDcol, masscol, xcol, ycol, zcol, Rvircol, bcol, ccol)
	print(len(LGID))
	
	omx, omy, omz = xyz_reader2(filename_high_res_filtered_overmassive, xcol, ycol, zcol)
#	#come up with fake binaries that are similar to the actual binaries from bringing two isolated haloes together
#	#produces a population file for the 1,858 pairs of fake binaries
	
#	truISOLG_finder()

	iso_candidate = loadtxt(isolated_lg_file_at_a_Dsep, dtype = int, usecols = (0,), unpack = True, skiprows = 1)
	x = np.where(iso_candidate == -1)
	print(len(x[0]))
	sys.exit()
	
	fakebin_finder()
	
	totm, totx, toty, totz = loadtxt(filename_high_res, skiprows = 1, usecols = (masscol, xcol, ycol, zcol), unpack=True)
	print("unpacked total halos")
	
	truISOLG_sat_finder()

	fakepair_pop()
	
	sys.exit()
	
	
#	iso_jlst, realhalo = loadtxt(Trueisolated_lg_file, unpack = True, skiprows = 1)
#	
#	readpopinfo2(isolated_lg_file_nosub, filename_trueISO_pop, iso_jlst)

#	#true isolated haloes info
#	
#	isoTLGID, isoTlgbm, isoTlgx, isoTlgy, isoTlgz, isoTlgRvir, isoTlgb, isoTlgc = col_reader5(filename_trueISO_pop, IDcol, masscol, xcol, ycol, zcol, Rvircol, bcol, ccol)
	
####################################################################################################################################
###	
######################## PLOT ISOLATED LG-MASS HALO POPULATION CHARACTERISTICS ####################################################
#	
##	print ("number of isolated lg mass haloes:", len(isoTlgbm))
###	histpopplot(isoTlgbm, "M [$h^{-1}$M$_{\odot}$]", "N", paperaddress("iso_mass.pdf"),  1, 0, 0, 20, "(a)")
##	print ("Histograph of LG mass distribution is saved as " + "iso_mass.pdf")
##	histpopplot(isoTlgRvir, "$R_{\matm{vir}} [\matm{h}^{-1}\matm{kpc}]$", "N", paperaddress("iso_virR.pdf"),  0, 0, 0, 40)
###	print "Histograph of LG virR distribution is saved as " + "iso_virR.pdf"
##	histpopplot_ecce(isoTlgb, isoTlgc, "axes ratio", "N", paperaddress("iso_bc.pdf"), 20, "(b)")
##	print ("Histograph of LG eccentricity b and c distribution is saved as " + "iso_bc.pdf")
##		
##	sys.exit()
######################### COMPOSE FAKE BINARY LG-MASS HALO  #########################################################################


