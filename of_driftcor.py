#! /usr/bin/python
# -*- coding: utf-8 -*-

#take the data interpolated to 12UTC (first column is still original time)
#interpolate spatially the ALS value
#interpolate time value

#extract for each virtual buoy position, time value
#calculate velocity for each virtual buoy


import numpy as np
from pyresample import image, geometry, kd_tree
import pyresample as pr
import matplotlib.pyplot as plt
from shapely.geometry import Polygon as shpol
from descartes import PolygonPatch
from matplotlib.path import Path
import csv
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
from datetime import date, timedelta, datetime
from of_func import *

path = '../data/frm_DTU_020816/ALS/'
#path = '../data/frm_DTU_040417/'
outpath_data = '../data/'

#ALS1
alsfile_of = path+'ALS_20150419_1032_5x5_corr_sh.frs'
#alsfile_of = path+'ALS_20150419_1032_5x5_raw_shr.frs'
date = datetime(2015,4,19,0)
date = datetime(2015,4,19,12)
gpsfile=outpath_data+'19042015nav.txt'
outname = path+'ALS_20150419_1032_5x5_corr_sh_dc.frs'
#outname = path+'ALS_20150419_1032_5x5_raw_dc.frs'

#ALS2
alsfile_of = path+'ALS_20150424_0955_5x5_corr_sh.frs'
#alsfile_of = path+'ALS_20150424_0955_5x5_raw.frs'
date = datetime(2015,4,24,0)
date2 = datetime(2015,4,24,12)
gpsfile=outpath_data+'24042015nav.txt'
outname = path+'ALS_20150424_0955_5x5_corr_sh_dc.frs'
#outname = path+'ALS_20150424_0955_5x5_raw_dc.frs'

#get gps data
f = open(gpsfile, 'r')
tmp = f.readlines()
#print help(f.readlines)
tmp_lst = []
for i in range(0,len(tmp)):
  line = tmp[i].replace('\r\n','').split(',')
  tmp_lst.extend(line)
  
#write gps file with structured rows and cols
n = 47
tmp_lst = np.array(tmp_lst).reshape(len(tmp_lst)/n,n)
#import pdb; pdb.set_trace()
#read all useful cols
gps_tm = np.array(map(lambda x: date+timedelta(seconds=np.float(x[:2])*3600+np.float(x[2:4])*60+np.float(x[4:])), tmp_lst[:,18]))	#UTC time on format hhmmss.ss
gps_lons = np.array(map(lambda x: np.float(x[0:3]) + np.float(x[3:])/60., tmp_lst[:,28]))
gps_lats = np.array(map(lambda x: np.float(x[0:2]) + np.float(x[2:])/60., tmp_lst[:,26]))

#read original als data
tmp = getColumn(alsfile_of,0,delimiter=' ')
als_tm = np.array(map(lambda x: date+timedelta(hours=float(x)),tmp))
als_fb = np.array(getColumn(alsfile_of,3,delimiter=' '),dtype=float)	#use column 3 for processed FB, clumn 6 for raw elevations
als_lats = np.array(getColumn(alsfile_of,1,delimiter=' '),dtype=float)
als_lons = np.array(getColumn(alsfile_of,2,delimiter=' '),dtype=float)
#als_elev = np.array(getColumn(alsfile_of,6,delimiter=' '),dtype=float)

#shorten the GPS arrays
start = als_tm[0]
end = als_tm[-1]
si = np.argmin(np.abs(gps_tm-start))
ei = np.argmin(np.abs(gps_tm-end))
print start
print end
print 'processing '+str(als_tm.shape[0])+' time steps'

gps_tm = gps_tm[si:ei]
gps_lons = gps_lons[si:ei]
gps_lats = gps_lats[si:ei]

#project positions
from pyproj import Proj,Geod
#g = Geod(ellps='WGS84')
epsg3575 = Proj("+init=EPSG:3575")

als_x,als_y = epsg3575(als_lons,als_lats)
gps_x,gps_y = epsg3575(gps_lons,gps_lats)

###look if data is spiky
##plt.plot(gps_tm,gps_lons)
#plt.plot(als_tm,als_lons)

##plt.plot(gps_tm,gps_y)
##plt.plot(als_tm,als_y)
#plt.show()
#exit()

idx12 = np.argmin(np.abs(gps_tm-date2))
print idx12
#calculate distance and direction of drift for every point in als file to 12UTC
als_x_dc = []
als_y_dc = []
x_shift = []
y_shift = []
for i in range(0,als_tm.shape[0]):
  idx = np.argmin(np.abs(gps_tm-als_tm[i]))
  xshift = gps_x[idx12] - gps_x[idx]
  yshift = gps_y[idx12] - gps_y[idx]
  print i, als_tm[i], xshift, yshift
  als_x_dc.append(als_x[i]+xshift)
  als_y_dc.append(als_y[i]+yshift)
  x_shift.append(xshift)
  y_shift.append(yshift)

als_x_dc = np.array(als_x_dc)
als_y_dc = np.array(als_y_dc)
x_shift = np.array(x_shift)
y_shift = np.array(y_shift)

#get also unprojected coordinates, just in case...
lonlat = Proj('+init=EPSG:4326')
from pyproj import transform
als_lons_dc,als_lats_dc = transform(epsg3575,lonlat,als_x_dc,als_y_dc)

#save a text file (table)
tt = [als_tm,als_lats,als_lats_dc,als_lons,als_lons_dc,als_y,als_x,als_y_dc,als_x_dc,y_shift,x_shift,als_fb]#,als_elev]
table = zip(*tt)

with open(outname, 'wb') as f:
  #header
  f.write(b'time, lat, dc lat, lon, dc lon, y, x, dc y, dc x, yshift (m), xshift (m), freeboard (m)\n')#, elevation (m)\n')
  np.savetxt(f, table, fmt="%s", delimiter=",")


