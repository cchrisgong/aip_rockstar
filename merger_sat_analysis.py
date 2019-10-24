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
import math
from numpy import amax, amin, median, sin, cos, pi, arccos

def RockstarLGB_MMAXPROG_path(Rockstarfilelist, filenum):
	filename_high_res = Rockstarfilelist[filenum]
	a = float(filename_high_res[6:12]) #scale
	z_redshift = (1.- a)/a #redshift
	# for each redshift there is a folder to store data in
	datfoldr = "LGB_merger/snapshot_for_z=%.4f" % z_redshift
	LGB_path = datfoldr + "/z_=%.4f_mmax.dat" % z_redshift
	return datfoldr, LGB_path, z_redshift

def sat_feature_calc(index_Asat_list, index_Bsat_list, argument_list, z_redshift):
	Acosangle_sat_1z = np.array([])
	Asinangle_sat_1z = np.array([])
	Av_relhost_sat_1z = np.array([])
	Av_relhost_cosangle_sat_1z = np.array([])
	Av_relcent_cosangle_sat_1z = np.array([])
	Apv_relhost_cosangle_sat_1z = np.array([])
	Apv_relcent_cosangle_sat_1z = np.array([])
	Adist_sat_1z = np.array([])
	Am_sat_1z = np.array([])
	AID_sat_1z = np.array([])
	AdescID_sat_1z = np.array([])

	Bcosangle_sat_1z = np.array([])
	Bsinangle_sat_1z = np.array([])
	Bv_relhost_sat_1z = np.array([])
	Bv_relhost_cosangle_sat_1z = np.array([])
	Bv_relcent_cosangle_sat_1z = np.array([])
	Bpv_relhost_cosangle_sat_1z = np.array([])
	Bpv_relcent_cosangle_sat_1z = np.array([])
	Bdist_sat_1z = np.array([])
	Bm_sat_1z = np.array([])
	BID_sat_1z = np.array([])
	BdescID_sat_1z = np.array([])

	for i in range(len(index_Asat_list)):
		index_Asat = index_Asat_list[i]
		index_Bsat = index_Bsat_list[i]
		try:
			[ID, descID, x, y, z, m, vx, vy, vz, Apos, Bpos, Avel, Bvel, ABvel, lineBA, dsep, dsep0] = argument_list[i]
		except IndexError:
			for l in range(len(argument_list[i])):
				argument_list[i][l] = np.array([argument_list[i][l]])
			[ID, descID, x, y, z, m, vx, vy, vz, Apos, Bpos, Avel, Bvel, ABvel, lineBA, dsep, dsep0] = argument_list[i]
		
		lineAB = -lineBA
		pos = np.transpose( np.array([x, y, z]) )
	
		Asatpos = np.transpose( np.array([x[index_Asat], y[index_Asat], z[index_Asat]]) ) 
		Bsatpos = np.transpose( np.array([x[index_Bsat], y[index_Bsat], z[index_Bsat]]) )
	
		AsattoA = coordsubstract_arr(Asatpos, Apos)
		BsattoB = coordsubstract_arr(Bsatpos, Bpos)
	
		AsattoA_norm = np.transpose(np.transpose(AsattoA)/np.linalg.norm(AsattoA, axis=1))
		BsattoB_norm = np.transpose(np.transpose(BsattoB)/np.linalg.norm(BsattoB, axis=1))
	
		Acosangle_sat_1z = np.append(Acosangle_sat_1z, np.dot( AsattoA_norm, lineAB) ) #angle of satellite position
		Bcosangle_sat_1z = np.append(Bcosangle_sat_1z, np.dot( BsattoB_norm, lineBA) )

		Asinangle_sat_1z = np.append(Asinangle_sat_1z, np.linalg.norm(np.cross(AsattoA_norm, lineAB), axis=1) ) #angle of satellite position
		Bsinangle_sat_1z = np.append(Bsinangle_sat_1z, np.linalg.norm(np.cross(BsattoB_norm, lineBA), axis=1) )
		
		Asatvel = np.transpose( np.array([vx[index_Asat], vy[index_Asat], vz[index_Asat]]))
		Bsatvel = np.transpose( np.array([vx[index_Bsat], vy[index_Bsat], vz[index_Bsat]]))
		
		Asatvel_reltohost = Asatvel - Avel
		Bsatvel_reltohost = Bsatvel - Bvel
		
		Asatvel_reltocent = Asatvel - ABvel
		Bsatvel_reltocent = Bsatvel - ABvel
		
		Av_relhost_sat_1z = np.append(Av_relhost_sat_1z, np.linalg.norm(Asatvel_reltohost, axis=1))
		Bv_relhost_sat_1z = np.append(Bv_relhost_sat_1z, np.linalg.norm(Bsatvel_reltohost, axis=1))
		
		Asatvel_norm = np.transpose(np.transpose(Asatvel)/np.linalg.norm(Asatvel, axis=1))
		Bsatvel_norm = np.transpose(np.transpose(Bsatvel)/np.linalg.norm(Bsatvel, axis=1))
		
		Asatvel_reltohost_norm = np.transpose(np.transpose(Asatvel_reltohost)/np.linalg.norm(Asatvel_reltohost, axis=1))
		Bsatvel_reltohost_norm = np.transpose(np.transpose(Bsatvel_reltohost)/np.linalg.norm(Bsatvel_reltohost, axis=1))
		
		Asatvel_reltocent_norm = np.transpose(np.transpose(Asatvel_reltocent)/np.linalg.norm(Asatvel_reltocent, axis=1))
		Bsatvel_reltocent_norm = np.transpose(np.transpose(Bsatvel_reltocent)/np.linalg.norm(Bsatvel_reltocent, axis=1))

		Av_relhost_cosangle_sat_1z = np.append(Av_relhost_cosangle_sat_1z, np.dot(Asatvel_reltohost_norm, lineAB) ) #angle of satellite velocity relative to host
		Bv_relhost_cosangle_sat_1z = np.append(Bv_relhost_cosangle_sat_1z, np.dot(Bsatvel_reltohost_norm, lineBA) )
		
		Av_relcent_cosangle_sat_1z = np.append(Av_relcent_cosangle_sat_1z, np.dot(Asatvel_reltocent_norm, lineAB) ) #angle of satellite velocity relative to center
		Bv_relcent_cosangle_sat_1z = np.append(Bv_relcent_cosangle_sat_1z, np.dot(Bsatvel_reltocent_norm, lineBA) )

		Apv_relhost_cosangle_sat_1z = np.append(Apv_relhost_cosangle_sat_1z, np.einsum('ij,ij->i', Asatvel_reltohost_norm, AsattoA_norm) ) #angle between position and velocity
		Bpv_relhost_cosangle_sat_1z = np.append(Bpv_relhost_cosangle_sat_1z, np.einsum('ij,ij->i', Bsatvel_reltohost_norm, BsattoB_norm) )

		Apv_relcent_cosangle_sat_1z = np.append(Apv_relcent_cosangle_sat_1z, np.einsum('ij,ij->i', Asatvel_reltocent_norm, AsattoA_norm) ) #angle between position and velocity
		Bpv_relcent_cosangle_sat_1z = np.append(Bpv_relcent_cosangle_sat_1z, np.einsum('ij,ij->i', Bsatvel_reltocent_norm, BsattoB_norm) )

		Adist_sat_1z = np.append( Adist_sat_1z, coorddistance_arr(pos[index_Asat], np.array([Apos]) )/dsep0 )
		Bdist_sat_1z = np.append( Bdist_sat_1z, coorddistance_arr(pos[index_Bsat], np.array([Bpos]) )/dsep0 )

		Am_sat_1z = np.append( Am_sat_1z, m[index_Asat] )
		Bm_sat_1z = np.append( Bm_sat_1z, m[index_Bsat] )

		AID_sat_1z = np.append( AID_sat_1z, ID[index_Asat] )
		BID_sat_1z = np.append( BID_sat_1z, ID[index_Bsat] )

		AdescID_sat_1z = np.append(AdescID_sat_1z, descID[index_Asat])
		BdescID_sat_1z = np.append(BdescID_sat_1z, descID[index_Bsat])
		
	Acosangle_sat_1z = Acosangle_sat_1z.flatten()
	Av_relhost_sat_1z = Av_relhost_sat_1z.flatten()
	Av_relhost_cosangle_sat_1z = Av_relhost_cosangle_sat_1z.flatten()
	Av_relcent_cosangle_sat_1z = Av_relcent_cosangle_sat_1z.flatten()
	Adist_sat_1z = Adist_sat_1z.flatten()
	Am_sat_1z = Am_sat_1z.flatten()
	Apv_relhost_cosangle_sat_1z = Apv_relhost_cosangle_sat_1z.flatten()
	Apv_relcent_cosangle_sat_1z = Apv_relcent_cosangle_sat_1z.flatten()
	AID_sat_1z = AID_sat_1z.flatten()
	AdescID_sat_1z = AdescID_sat_1z.flatten()
	Asinangle_sat_1z = Asinangle_sat_1z.flatten()
	
	Bcosangle_sat_1z = Bcosangle_sat_1z.flatten()
	Bv_relhost_sat_1z = Bv_relhost_sat_1z.flatten()
	Bv_relhost_cosangle_sat_1z = Bv_relhost_cosangle_sat_1z.flatten()
	Bv_relcent_cosangle_sat_1z = Bv_relcent_cosangle_sat_1z.flatten()
	Bdist_sat_1z = Bdist_sat_1z.flatten()
	Bm_sat_1z = Bm_sat_1z.flatten()
	Bpv_relhost_cosangle_sat_1z = Bpv_relhost_cosangle_sat_1z.flatten()
	Bpv_relcent_cosangle_sat_1z = Bpv_relcent_cosangle_sat_1z.flatten()
	BID_sat_1z = BID_sat_1z.flatten()
	BdescID_sat_1z = BdescID_sat_1z.flatten()
	Bsinangle_sat_1z = Bsinangle_sat_1z.flatten()

	DAT = np.column_stack((AID_sat_1z, AdescID_sat_1z, Acosangle_sat_1z, Av_relhost_cosangle_sat_1z, Av_relcent_cosangle_sat_1z, Apv_relhost_cosangle_sat_1z, Apv_relcent_cosangle_sat_1z, Adist_sat_1z, Am_sat_1z, Av_relhost_sat_1z, Asinangle_sat_1z))
	np.savetxt('z_=%.4f_sat_merger_samesizesample_hostA.dat' % z_redshift, DAT, delimiter=" ")
	DAT = np.column_stack((BID_sat_1z, BdescID_sat_1z, Bcosangle_sat_1z, Bv_relhost_cosangle_sat_1z, Bv_relcent_cosangle_sat_1z, Bpv_relhost_cosangle_sat_1z, Bpv_relcent_cosangle_sat_1z, Bdist_sat_1z, Bm_sat_1z, Bv_relhost_sat_1z, Bsinangle_sat_1z))
	np.savetxt('z_=%.4f_sat_merger_samesizesample_hostB.dat' % z_redshift, DAT, delimiter=" ")
			
def merger_sat_path(filenum):
	filename_high_res = Rockstarfilelist[filenum]
#	print(filename_high_res)
	a = float(filename_high_res[6:12]) #scale
	z_redshift = (1.- a)/a #redshift
	datmainfoldr = "/LGB_merger"
	datfoldr = "/snapshot_for_z=%.4f" % z_redshift
	
	if filenum == 0:
		satfoldrz0 = "/store/erebos/cgong/Rockstar_smaller_massrange" + datmainfoldr + datfoldr + "/smallhalos_highres" #folder of z = 0 sat files
		sat_files = os.listdir(satfoldrz0)
		sat_files.sort(key = natural_keys)
		datfoldrz0 = "snapshot%.4f" % a
		neglist = loadtxt(pwd0 + "/" + datfoldrz0 + "/negfile.dat")
		return satfoldrz0, neglist, sat_files, z_redshift
	
	else:
		# merger tree satellites
		LGBsat_path = "/store/erebos/cgong/Rockstar_smaller_massrange" + datmainfoldr + datfoldr + "/smallhalos_highres_merg"
		sat_files = os.listdir(LGBsat_path)
		sat_files.sort(key = natural_keys)
		
		return LGBsat_path, sat_files, z_redshift
	
def merger_sat_analysis_outputID():
	maxmass_ID_Asat = []
	maxmass_ID_Bsat = []
	maxmass_descID_Asat = []
	maxmass_descID_Bsat = []
	maxmass_mass_Asat = []
	maxmass_mass_Bsat = []	
	##########################################################################################################
	# z=0 zero redshift satellite: 
	##########################################################################################################
	# host
	datfoldr, hostfile0, z0 = RockstarLGB_MMAXPROG_path(Rockstarfilelist, 0)
	IDhost, mhost, vxhost, vyhost, vzhost, xhost, yhost, zhost, rvirhost, bhost, chost = loadtxt(hostfile0, unpack =1)
	massivehost = []
	smallerhost = []
	for i in range(len(IDhost)): # sort hosts according to mass
		if i % 2 == 0:
			hostAmass = mhost[i]
			hostBmass = mhost[i + 1]
			if hostAmass > hostBmass:
				massivehost.append(i)
				smallerhost.append(i + 1)
			else:
				massivehost.append(i + 1)
				smallerhost.append(i)
	
	hostApos = np.transpose(np.array([xhost[massivehost], yhost[massivehost], zhost[massivehost]]))
	hostBpos = np.transpose(np.array([xhost[smallerhost], yhost[smallerhost], zhost[smallerhost]]))
	dsep = coorddistance_arr(hostApos, hostBpos)
	
	hostAID = IDhost[massivehost]
	hostBID = IDhost[smallerhost]
	
	# satellite		
	satfoldrz0, neglist, sat_files, z_redshift = merger_sat_path(0)
	IDold_satA_1z = []
	IDold_satB_1z = []
	descIDold_satA_1z = []
	descIDold_satB_1z = []	
	mold_satA_1z = []
	mold_satB_1z = []

	for i in range(int(len(xhost)/2)):
		if i not in neglist:
			
			IDold, descIDold, mold, xold, yold, zold = col_reader8( satfoldrz0 + '/' + str(i) + '_small_halos_feats.dat', IDcol, descIDcol, masscol, xcol, ycol, zcol ) #sat of last redshift
			pos = np.transpose( np.array([xold, yold, zold]) )
			
			rhalf = dsep[i] / 2e0
			
			Apos = hostApos[i]
			Bpos = hostBpos[i]
			
			index_Asat = np.where( (coorddistance_arr(pos, np.array([Apos])) < rhalf) & ( IDold!= hostAID[i] ))
			index_Bsat = np.where( (coorddistance_arr(pos, np.array([Bpos])) < rhalf) & ( IDold!= hostBID[i] ))

			IDold_satA_1z.append( IDold[index_Asat[0]] ) #all the IDs of A's satellites one redshift earlier
			IDold_satB_1z.append( IDold[index_Bsat[0]] )
			descIDold_satA_1z.append( descIDold[index_Asat[0]] )
			descIDold_satB_1z.append( descIDold[index_Bsat[0]] )
			mold_satA_1z.append(mold[index_Asat[0]])
			mold_satB_1z.append(mold[index_Bsat[0]])
			
	maxmass_ID_Asat.append(IDold_satA_1z)
	maxmass_ID_Bsat.append(IDold_satB_1z)
	maxmass_mass_Asat.append(mold_satA_1z)
	maxmass_mass_Bsat.append(mold_satB_1z)
	maxmass_descID_Asat.append(descIDold_satA_1z)
	maxmass_descID_Bsat.append(descIDold_satB_1z)
	
	##########################################################################################################
	# z>0 zero redshift satellite:
	##########################################################################################################
	# host
	
	for l in range(1, 27):
#	for l in range(1, 3):
	#		print (l)
		datfoldr, hostfile, z_redshift = RockstarLGB_MMAXPROG_path(Rockstarfilelist, l)
		LGBsat_path, sat_files, z_redshift = merger_sat_path(l)

		IDhost, mhost, vxhost, vyhost, vzhost, xhost, yhost, zhost, rvirhost, bhost, chost = loadtxt(hostfile, unpack =1)
		
		IDnew_satA_1z = []
		descIDnew_satA_1z =[]
		mnew_satA_1z = []
		mnew_satB_1z = []
		IDnew_satB_1z = []
		descIDnew_satB_1z =[]
		
		counter = 0 # represents one host pair ( needs to be counted because some are in negative list (due to large satellite))
	
		for i in range(int(len(xhost)/2)):
#		for i in range(100):
			if i not in neglist:
#				print(str(i) + '_small_halos_feats.dat')
				IDnew, descIDnew, m, x, y, z, vx, vy, vz = col_reader7(LGBsat_path + '/' + str(i) + '_small_halos_feats.dat', IDcol, descIDcol, masscol, xcol, ycol, zcol, vxcol, vycol, vzcol)
				
				index_Asat = []
				index_Bsat = []
				
#				cas = np.where(descIDnew == -1) #no Descendent # this is not needed, as the code in merger_sat.py garantees that all trace backs have a descendent in present day z=0
				
				for ID_oldA in maxmass_ID_Asat[-1][counter]: #for each pair's one single satellite at the last time slice
					index_Asat0 = np.where( descIDnew == ID_oldA ) #All A's satellite's progenitors
					index_Asat0 = index_Asat0[0]
					if len(index_Asat0) != 0:
						try:
							mmax = np.amax(m[index_Asat0]) #pick out biggest predecessor mass
							j = np.where(m[index_Asat0] == mmax)
							index_Asat.append(int(index_Asat0[j[0]][0]))
						except IndexError: #if index_Asat0 is a scalar, i.e. only one progenitor present
							mmax = m
							index_Asat.append(0)
							
					else:  #if cannot find progenitor in the last time slice, No progenitor
#						if ID_oldA ==810807520 or ID_oldA == 810806275:
						for n in range(len(maxmass_ID_Asat))[::-1]:
							bas = np.where(maxmass_ID_Asat[n][counter] == ID_oldA)
							ID_oldA = maxmass_descID_Asat[n][counter][bas[0]]
							maxmass_ID_Asat[n][counter] = np.delete(maxmass_ID_Asat[n][counter], bas[0])
							maxmass_descID_Asat[n][counter] = np.delete(maxmass_descID_Asat[n][counter], bas[0])
							maxmass_mass_Asat[n][counter] = np.delete(maxmass_mass_Asat[n][counter], bas[0])
							
				for ID_oldB in maxmass_ID_Bsat[-1][counter]: #for each pair's one single satellite at the last time slice
					index_Bsat0 = np.where( descIDnew == ID_oldB ) #All B's satellite's progenitors
					index_Bsat0 = index_Bsat0[0]
					if len(index_Bsat0) != 0:
						try:
							mmax = np.amax(m[index_Bsat0]) #pick out biggest predecessor mass
							j = np.where(m[index_Bsat0] == mmax)
							index_Bsat.append(int(index_Bsat0[j[0]][0]))
						except IndexError:
							mmax = m
							index_Bsat.append(0)
							
					else: # if cannot find progenitor in the last time slice, i.e. running into a dead end, we need to kick the descendent out 
						for n in range(len(maxmass_ID_Bsat))[::-1]:
							bas = np.where(maxmass_ID_Bsat[n][counter] == ID_oldB)
							ID_oldB = maxmass_descID_Bsat[n][counter][bas[0]]
							maxmass_ID_Bsat[n][counter] = np.delete(maxmass_ID_Bsat[n][counter], bas[0])							
							maxmass_descID_Bsat[n][counter] = np.delete(maxmass_descID_Bsat[n][counter], bas[0])
							maxmass_mass_Bsat[n][counter] = np.delete(maxmass_mass_Bsat[n][counter], bas[0])				
				
				counter += 1 #because there are some binaries on the neglist due to their large satellite (>.5 minor host mass)
				
				try:
					IDnew_satA_1z.append(IDnew[index_Asat])
					descIDnew_satA_1z.append(descIDnew[index_Asat])
					mnew_satA_1z.append(m[index_Asat])
					
				except IndexError:
					IDnew_satA_1z.append([IDnew])
					descIDnew_satA_1z.append([descIDnew])
					mnew_satA_1z.append([m])
				try:
					IDnew_satB_1z.append(IDnew[index_Bsat])
					descIDnew_satB_1z.append(descIDnew[index_Bsat])
					mnew_satB_1z.append(m[index_Bsat])
					
				except IndexError:
					IDnew_satB_1z.append([IDnew])
					descIDnew_satB_1z.append([descIDnew])
					mnew_satB_1z.append([m])

		maxmass_descID_Asat.append(descIDnew_satA_1z)
		maxmass_descID_Bsat.append(descIDnew_satB_1z)
		maxmass_ID_Asat.append(IDnew_satA_1z)
		maxmass_ID_Bsat.append(IDnew_satB_1z)
		maxmass_mass_Asat.append(mnew_satA_1z)
		maxmass_mass_Bsat.append(mnew_satB_1z)				
		IDold_satA_1z = IDnew_satA_1z
		IDold_satB_1z = IDnew_satB_1z
#		for n in range(len(maxmass_ID_Asat))[::-1]:
#			print(n, len(maxmass_ID_Asat[n][0]))
		
	for k in range(len(maxmass_ID_Asat)):
		counter = 0
		datfoldr, hostfile, z_redshift = RockstarLGB_MMAXPROG_path(Rockstarfilelist, k)
		os.chdir(datfoldr)
		datfoldr_sat = "sat_analysis"
		try:
			os.makedirs(datfoldr_sat)
		except OSError:
			pass		
		for i in range(int(len(xhost)/2)):
			if i not in neglist:
#				print(len(maxmass_ID_Asat[k][counter]))
				np.savetxt('sat_analysis/z=%.4f_satID_hostA%d.dat' % (z_redshift,i), maxmass_ID_Asat[k][counter], fmt='%i', delimiter=" ")
				np.savetxt('sat_analysis/z=%.4f_satID_hostB%d.dat' % (z_redshift,i), maxmass_ID_Bsat[k][counter], fmt='%i', delimiter=" ")			
				counter += 1
		os.chdir("../..")

def merger_sat_analysis():
	massivehost = []
	smallerhost = []
	
	for l in range(27):
		#sat files path
		if l ==0:
			satfoldr, neglist, sat_files, z_redshift = merger_sat_path(l)
		else:
			satfoldr, sat_files, z_redshift = merger_sat_path(l)
		
		#host file path
		datfoldr, hostfile, z = RockstarLGB_MMAXPROG_path(Rockstarfilelist, l)
		IDhost, mhost, vxhost, vyhost, vzhost, xhost, yhost, zhost, rvirhost, bhost, chost = loadtxt(hostfile, unpack =1)

		if l == 0:
			for i in range(len(IDhost)):
				if i %2 == 0:
					if int(i/2) not in neglist:
						hostAmass = mhost[i]
						hostBmass = mhost[i + 1]
#						print('A', xhost[i],yhost[i],zhost[i])
#						print('B',xhost[i+1],yhost[i+1],zhost[i+1])
						if hostAmass > hostBmass:
							massivehost.append(i)
							smallerhost.append(i + 1)
						else:
							massivehost.append(i + 1)
							smallerhost.append(i)
			
			hostAID = np.transpose(np.array([IDhost[massivehost], IDhost[massivehost], IDhost[massivehost]]))
			hostBID = np.transpose(np.array([IDhost[smallerhost], IDhost[smallerhost], IDhost[smallerhost]]))

		hostApos = np.transpose(np.array([xhost[massivehost], yhost[massivehost], zhost[massivehost]]))
		hostBpos = np.transpose(np.array([xhost[smallerhost], yhost[smallerhost], zhost[smallerhost]]))
		if l == 0:
			dsep0 = coorddistance_arr(hostApos, hostBpos)
		dsep = coorddistance_arr(hostApos, hostBpos)
		hostAvel = np.transpose(np.array([vxhost[massivehost], vyhost[massivehost], vzhost[massivehost]]))
		hostBvel = np.transpose(np.array([vxhost[smallerhost], vyhost[smallerhost], vzhost[smallerhost]]))
		hostABvel = (hostAvel + hostBvel) / 2e0
		
		argument_list = []
		index_Asat_list = []
		index_Bsat_list = []
		
		os.chdir(datfoldr)
		counter = 0
		
		for i in range(int(len(xhost)/2)):
			if i not in neglist:
#				print(str(i) + '_small_halos_feats.dat')
				IDnew, descIDnew, m, x, y, z, vx, vy, vz = col_reader7( satfoldr + '/' + str(i) + '_small_halos_feats.dat', IDcol, descIDcol, masscol, xcol, ycol, zcol, vxcol, vycol, vzcol ) #sat of last redshift
				pos = np.transpose( np.array([x, y, z]) )
				
				Apos = hostApos[counter]
				Bpos = hostBpos[counter]
				lineBA = coordsubstract(Apos, Bpos)/dsep[counter] # vector B pointing towards A
				
				Avel = hostAvel[counter]
				Bvel = hostBvel[counter]
				ABvel = hostABvel[counter]
				
				argument = [IDnew, descIDnew, x, y, z, m, vx, vy, vz, Apos, Bpos, Avel, Bvel, ABvel, lineBA, dsep[counter], dsep0[counter]]
				argument_list.append(argument)
				
				ID_satA_1z = loadtxt('sat_analysis/z=%.4f_satID_hostA%d.dat' % (z_redshift, i), dtype = 'int') 
				ID_satB_1z = loadtxt('sat_analysis/z=%.4f_satID_hostB%d.dat' % (z_redshift, i), dtype = 'int')
		
				index_Asat = np.where( np.in1d( IDnew, ID_satA_1z ) )
				index_Bsat = np.where( np.in1d( IDnew, ID_satB_1z ) )
				index_Asat_list.append(np.array(index_Asat[0]))
				index_Bsat_list.append(np.array(index_Bsat[0]))
				counter += 1
		sat_feature_calc(index_Asat_list, index_Bsat_list, argument_list, z_redshift)
		os.chdir("../..")
						
if __name__ == "__main__":
#	 gathering all the catalog files into a list, ordered by snapshots from current day snapshots to beginning of the universe,

	Rockstarfilelist0 = []
	Rockstar_files = os.listdir(Rockstar_hlist_address)
	pwd0 = os.getcwd()

	for i in range(len(Rockstar_files)):
		filename = Rockstar_files[i]
		if filename.endswith("list"):
			Rockstarfilelist0.append(filename)
	Rockstarfilelist0.sort(key = natural_keys)
	Rockstarfilelist = Rockstarfilelist0[:][::-1]
	merger_sat_analysis_outputID()
	merger_sat_analysis()
