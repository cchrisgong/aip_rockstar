'''
non-projected: finding small halos by locating them in the 3D box spanned by dsep surrounding the midpoint of a LG binary
projected: finding small halos by locating them in the 2D square spanned by dsep surrounding the midpoint of a LG binary
'''

import os
from helpers import *
from constant import *
from numpy import loadtxt
import numpy as np
import sys

def nonindexed_smallhaloFinder(filename_high_res, LGID, LGm, LGx, LGy, LGz, SMALLm, SMALLx, SMALLy, SMALLz):
	#This function takes all the list of binaries without subs and search for small satellites near the pair
	#This also throws out the pair when satellites appear too big, >.5 minor host
	#SMALLm, SMALLx, SMALLy, SMALLz, SMALLID includes all haloes, not just the small ones according to range definition 
	f = open(filename_high_res,"r") #read header
	lines = f.readlines()
	f.close()
	
	f2 = open("negfile.dat","w")
	
	for j in range(len(LGID)): #for each LG halo
		
		if j % 2 == 0: #for every other LG halo (i.e. for every pair of LG halo)
			pairno = int(j / 2e0)
#			if pairno > 9919:
#			print ("small halo finder for pair", pairno)
			small_halo_indices_filename = str(pairno) + "_small_halos_feats.dat"
			small_halo_path = os.path.join(os.getcwd(), "smallhalos_highres", small_halo_indices_filename) 
#			if os.path.isfile(small_halo_path): #if file exists
#				pass
#			else:	
			f = open(small_halo_path, "w")
			hostApos = np.array([LGx[j], LGy[j], LGz[j]])
			hostBpos = np.array([LGx[j+1], LGy[j+1], LGz[j+1]])

			hostAID = LGID[j]
			hostBID = LGID[j+1]

			midpos = np.array(middlecoordmpc(hostApos, hostBpos))
			dsep = coorddistance(hostApos, hostBpos)
			minorhost = min([LGm[j], LGm[j+1]])

			#construct 3 x 3 x 3 box around the midpoint of the binary with wrapped coordinates so they fall into simulation box
			boxmin = coordwrap2(midpos - dsep) #within a box extended by 1.5Mpc * 2 in three direction with the jth halo 
			boxmax = coordwrap2(midpos + dsep) #at its center
			SMALL = np.transpose(np.array([SMALLx, SMALLy, SMALLz]))

			f.write(lines[0])
			index_sat_inbox = np.where( np.ndarray.all(np.mod(SMALL-boxmin-dim/2e0, dim)-dim/2e0 > 0, axis = 1) & np.ndarray.all(np.mod(SMALL-boxmax-dim/2e0, dim)-dim/2e0 < 0, axis = 1) )
			index_sat = np.where( (coorddistance_arr(SMALL[index_sat_inbox[0]], np.array([midpos])) < dsep) & (coorddistance_arr(SMALL[index_sat_inbox[0]], hostApos) != 0 ) & (coorddistance_arr(SMALL[index_sat_inbox[0]], hostBpos) != 0 ) )
#			print(SMALL[index_Asat_inbox[0]][index_Asat[0]], coorddistance_arr(SMALL[index_Asat_inbox[0]], np.array([midpos])))
			for i in index_sat_inbox[0][index_sat[0]]:
				if SMALLm[i] > .5 * minorhost and SMALLm[i] < minorhost:
					f2.write(str(pairno) + "\n")
					f.close()
					os.remove(small_halo_path)
					break
				elif SMALLm[i] > minorhost:
					print (SMALLm[i], minorhost, SMALL[i], hostApos, hostBpos )
#					print ("something wrong!")
					f2.write(str(pairno) + "\n")
					f.close()
					os.remove(small_halo_path)					
#					sys.exit()
					break
				else:
					f.write(lines[i+headerlength])
			f.close()
		
	f2.close()

