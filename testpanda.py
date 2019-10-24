from numpy import loadtxt
import pandas as pd

import time
import os
from constant import *
from helpers import *

IDcol = 1
Rockstarfilelist0 = []
Rockstar_files = os.listdir(Rockstar_hlist_address)
pwd0 = os.getcwd()

for i in range(len(Rockstar_files)):
	filename = Rockstar_files[i]
	if filename.endswith("list"):
		Rockstarfilelist0.append(filename)
Rockstarfilelist0.sort(key = natural_keys)

Rockstarfilelist = Rockstarfilelist0[:][::-1]
#copyfile(filepath_catalog, newfilepath)
z0_filename = Rockstarfilelist[1] #z = 0 catalog file
print(z0_filename)
snapnum = z0_filename[6:12]
datfoldrz0 = "/snapshot" + snapnum 
z0_LGB_path = pwd0 + datfoldrz0 + "/" + z0_filename #z = 0 LGB file
start_time = time.time()
#ID = pd.read_csv(z0_LGB_path, skiprows = 16, sep = '\t', usecols = (IDcol,), dtype='i8') #loading IDs of binaries at z= 0
#ID = pd.read_csv(z0_LGB_path, sep = '\t', usecols = (1,), skiprows = headerlength) #loading IDs of binaries at z= 0
ID = pd.read_csv(z0_LGB_path, skiprows = headerlength, sep = ' ', usecols = (IDcol, descIDcol), dtype='i8')
time1 = time.time()
print(time1-start_time)
print (ID.head(), ID.head())

ID = loadtxt(z0_LGB_path, skiprows = 1, usecols = (IDcol,), unpack=True, dtype='i8') #loading IDs of binaries at z= 0
time2 = time.time()
print(time2-time1)
