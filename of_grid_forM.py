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
from of_func import *

path = '../data/frm_DTU_020816/ALS/'
outpath = '../plots/'

#ALS1
alsfile = path+'ALS_20150419_1032_5x5_lance_corr_12UTC.out'
vbuoys = '../vbuoys_1904_corr_added.csv'
reg = 'als1_zoom'
name_als = 'als19_fb_zoom.png'
name_time = 'als19_tm_zoom.png'

#ALS2
alsfile = path+'ALS_20150424_0955_5x5_lance_corr_12UTC.out'
#open water polygon
ow2 = '../open_water_24_04'
vbuoys = '../vbuoys_2404_corr_added.csv'
reg = 'als2_zoom'
name_als = 'als24_fb_zoom.png'
name_time = 'als24_tm_zoom.png'


#read data
als = np.array(getColumn(alsfile,3,delimiter=' '),dtype=float)
time = np.array(getColumn(alsfile,0,delimiter=' '),dtype=float)
lons = np.array(getColumn(alsfile,2,delimiter=' '),dtype=float)
lats = np.array(getColumn(alsfile,1,delimiter=' '),dtype=float)

np.save('lat_forM',lats)
np.save('lon_forM',lons)
exit()

#get the virtual buoy data
data = np.loadtxt(vbuoys, dtype=np.str,delimiter=',', skiprows=1)
latitude1 = np.ma.array(data[:,1],mask=data[:,1]=='0',dtype=float)
longitude1 = np.ma.array(data[:,0],mask=data[:,2]=='0',dtype=float)

#define projection and interpolate to regular grid
area_def = pr.utils.load_area('area.cfg', reg)
swath_def = geometry.SwathDefinition(lons=lons, lats=lats)
#als_map = kd_tree.resample_gauss(swath_def, als,area_def, radius_of_influence=10, sigmas=5)
als_map = kd_tree.resample_nearest(swath_def, als,area_def, radius_of_influence=10)
#time_map = kd_tree.resample_gauss(swath_def, time,area_def, radius_of_influence=10, sigmas=5)
time_map = kd_tree.resample_nearest(swath_def, time,area_def, radius_of_influence=10)

##mask irrelevant and extreme data
#mask = (time_map<10)|(time_map>14)|(als_map<-1)|(als_map>2)
#als_map = np.ma.array(als_map, mask=mask)
#time_map = np.ma.array(time_map, mask=mask)

##plotting
#pr.plot.save_quicklook(outpath+'als14_quick.png', area_def, result, label='ALS freeboard (m)')
#ALS data
bmap = pr.plot.area_def2basemap(area_def)
#bmng = bmap.bluemarble()
col = bmap.imshow(als_map, origin='upper')
x1, y1 = bmap(longitude1, latitude1)
bmap.plot(x1,y1,'o',linewidth=2,color='purple')
bmap.colorbar()
plt.savefig(outpath+name_als, bbox_inches='tight')

#time
bmap = pr.plot.area_def2basemap(area_def)
col = bmap.imshow(time_map, origin='upper')
x1, y1 = bmap(longitude1, latitude1)
bmap.plot(x1,y1,'o',linewidth=2,color='purple')
bmap.colorbar()
plt.savefig(outpath+name_time, bbox_inches='tight')

#get data at the buoy points
#print dir(area_def)
#print help(area_def.get_xy_from_lonlat)

#this is how i get indexes of buoy coordinates in the array of projected coordinates (x,y)
#why are some of the values masked? All buoys should be inside of the defined area...
#ok, now i know that positions of the first buoy are at least 777 m off in latitude. The rest might be similar distorted.
bix,biy = area_def.get_xy_from_lonlat(longitude1,latitude1)
print bix
print biy
print bix.shape
print bix[0],biy[0]
#this is how i get the whole array of projected area (gridded 1000 by 1000 points, 5m spacing)
x,y = area_def.get_proj_coords()
print x.shape
#this gives me projected x coordinate of the first buoy
print x[bix[0],biy[0]]
#how can i get a list (vector) of x,y coordinates for all buoys
#find nearest als value and time value for all the buoys (should have same index!)
bx = []
by = []
bals = []
btm = []
for i in range(0,68):
  print bix[i],biy[i]
  bx.append(x[bix[i],biy[i]])
  by.append(y[bix[i],biy[i]])
  bals.append(als_map[bix[i],biy[i]])
  btm.append(time_map[bix[i],biy[i]])

#double check where bx and by are on the map!
print bx
print by
#print bals
#print btm

print dir(area_def)
print help(area_def.get_lonlat)
print area_def.get_lonlat(bix[0],biy[0])
#outputs: (array(15.265880850705882), array(82.73487172127551))
print longitude1[0]
#outputs: 15.32202059
print latitude1[0]
#outputs:82.74126076

##700 m difference in latitude!!!
#(array(15.265880850705882), array(82.73487172127551))
#15.32202059
#82.74126076


##calculate displacements and velocities between 2 overflights


