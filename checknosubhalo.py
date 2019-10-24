from constant import *
import sys
def checknosubhalo_pair(tobefilteredfile, resultfile):
	f = open(tobefilteredfile, "r")
	lines = f.readlines()
	f.close()
	f1 = open(resultfile, "w")
	counter_num_subhalos = 0
	counter_num_notsubhalos = 0
	f1.write(lines[0])
	for i in range(1, len(lines)):
		if i % 2 == 0:
			s1 = lines[i-1].strip()
			datlist1 = [float(x) for x in s1.split()]
			s2 = lines[i].strip()
			datlist2 = [float(x) for x in s2.split()]
			if datlist1[1] == 0 and datlist2[1] == 0: #subhalo ID is not 0, i.e. it is a subhalo
				f1.write(lines[i-1])
				f1.write(lines[i])
				counter_num_notsubhalos += 1
			else:
				counter_num_subhalos += 1
				
	f1.close()
	return counter_num_notsubhalos
#	print (len(lines)-1)/2, counter_num_subhalos, counter_num_notsubhalos

def checknosubhalo(tobefilteredfile, resultfile, negresultfile):
	#write individual haloes into resultfile, write subhaloes into negresultfile
	subhaloID = loadtxt(tobefilteredfile, skiprows = 1, usecols = (subhalocol,), unpack=True)	
	
	counter_num_subhalos = 0
	f1 = open(resultfile, "w")
	f2 = open(negresultfile, "w")
	
	f = open(tobefilteredfile,"r") #read header
	line = f.readline()
	f.close()
	
	f1.write(line)# skip header
	f2.write(line)# skip header

	f = open(tobefilteredfile,"r") #read header
	lines = f.readlines()
	f.close()
	
	nonSHindex = np.where( subhaloID == -1 )
	for i in nonSHindex[0]:
		f1.write(lines[i+1])
	f1.close()
	
	SHindex = np.where( subhaloID != -1 )
	for i in SHindex[0]:
		f2.write(lines[i+1])
	f2.close()	
	
	return len(nonSHindex[0]), len(SHindex[0])
	
if __name__ == "__main__":
	checknosubhalo(filename_high_res_lg_wsub, filename_high_res_lg)
