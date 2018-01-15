#! /usr/bin/python
# -*- coding: utf-8 -*-

#export virtual buoys vector layer created in QGIS to CSV files with geographical projection (WGS84, EPSG:4326)
#plot the ALS data and points

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
import osr, gdal, ogr
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon as shpol
from shapely.ops import cascaded_union
from descartes import PolygonPatch
import os

path = '../data/'
outpath = '../plots/'

#als data
alsfile1 = '../data/ALS_20150419_EPSG3575.tiff'
alsfile2 = '../data/ALS_20150424_EPSG3575.tiff'
#virtual buoys
vbuoys1 = '../vbuoys_1904_corr_added.csv'
vbuoys2 = '../vbuoys_2404_corr_added.csv'

#######################################33
#1st overflight (19/04/2015)
outals1 = '../ALS_20150419'
outname1 = 'map_19_04_new'
outfile1 = 'area_19_04.shp'
outfile_tri1 = 'tri_19_04.shp'
title1 = 'ICE-ARC overflight 19/04/2015'
llcrnrlon1=13.8; llcrnrlat1=83.1
urcrnrlon1=14.25; urcrnrlat1=83.14

#2nd overflight (24/04/2015)
outals2 = '../ALS_20150424'
outname2 = 'map_24_04_new'
outfile2 = 'area_24_04.shp'
outfile_tri2 = 'tri_24_04.shp'
title2 = 'ICE-ARC overflight 24/04/2015'
llcrnrlon2=15.1; llcrnrlat2=82.72
urcrnrlon2=15.51; urcrnrlat2=82.761
########################################33

#get the virtual buoy data
data = np.loadtxt(vbuoys1, dtype=np.str,delimiter=',', skiprows=1)
latitude1 = np.ma.array(data[:,1],mask=data[:,1]=='0',dtype=float)
longitude1 = np.ma.array(data[:,0],mask=data[:,2]=='0',dtype=float)
#idb1 = np.array(data[:,2],dtype=int)

data = np.loadtxt(vbuoys2, dtype=np.str,delimiter=',', skiprows=1)
latitude2 = np.ma.array(data[:,1],mask=data[:,1]=='0',dtype=float)
longitude2 = np.ma.array(data[:,0],mask=data[:,2]=='0',dtype=float)
#idb2 = np.array(data[:,2],dtype=int)

m1 = Basemap(resolution='h', epsg = 3575,
	    #ellps = 'WGS84', projection='laea', lat_0=90., lon_0=10.,
	    llcrnrlon=llcrnrlon1, llcrnrlat=llcrnrlat1,
	    urcrnrlon=urcrnrlon1, urcrnrlat=urcrnrlat1)


m2 = Basemap(resolution='h', epsg = 3575,
	    #ellps = 'WGS84', projection='laea', lat_0=90., lon_0=10.,
	    llcrnrlon=llcrnrlon2, llcrnrlat=llcrnrlat2,
	    urcrnrlon=urcrnrlon2, urcrnrlat=urcrnrlat2)


#get the virtual buoys in geo projection
x1, y1 = m1(longitude1, latitude1)
x2, y2 = m2(longitude2, latitude2)

#get the ALS data
mv =-999

##first overflight
ds = gdal.Open(alsfile1)
als1 = ds.ReadAsArray()
#this will change all nan to mv and mask them at the same time
als1 = np.ma.fix_invalid(als1, fill_value=mv)
als1.mask = np.ma.nomask
#coords
xx1 = np.load('xx1.npy')
yy1 = np.load('yy1.npy')

##second overflight
ds = gdal.Open(alsfile2)
als2 = ds.ReadAsArray()
#this will change all nan to mv and mask them at the same time
als2 = np.ma.fix_invalid(als2, fill_value=mv)
als2.mask = np.ma.nomask
#coords
xx2 = np.load('xx2.npy')
yy2 = np.load('yy2.npy')

#get als values for each buoy

#als_vb1 = []
#als_vb2 = []
#for i in range(0,len(x1)):
  #print 'virtual buoy: ',i
  ##first overflight
  #box = ~((yy1 < y1[i]+3) & (yy1 > y1[i]-3) & (xx1 < x1[i]+3) & (xx1 > x1[i]-3) & (als1.T!=mv))
  #tmp = np.ma.array(als1.T,mask=box).compressed()
  #als_vb1.append(np.mean(tmp))
  ##second overflight
  #box = ~((yy2 < y2[i]+3) & (yy2 > y2[i]-3) & (xx2 < x2[i]+3) & (xx2 > x2[i]-3) & (als2.T!=mv))
  #tmp = np.ma.array(als2.T,mask=box).compressed()
  #als_vb2.append(np.mean(tmp))
  
#print  als_vb1 
#print  als_vb2

#np.save('als_vb1',als_vb1)
#np.save('als_vb2',als_vb2)

als_vb1 = np.load('als_vb1.npy')
als_vb2 = np.load('als_vb2.npy')


##calculate differences between overflights
#diff = np.array(als_vb2) - np.array(als_vb1)

##diff_thin = np.ma.array(diff,mask=als_vb1>0.2).compressed()
#diff_thick = np.ma.array(diff,mask=als_vb1<0.1).compressed()

##print diff
##print diff_thin

##plot a histogram of differences
#n, bins, patches = plt.hist(diff, 20, normed=0, facecolor='green', alpha=0.5, label='all')
##n, bins, patches = plt.hist(diff_thin, 20, normed=0, facecolor='red', alpha=0.5, label='thin')
#n, bins, patches = plt.hist(diff_thick, 20, normed=0, facecolor='blue', alpha=0.5, label='thick')

### add a 'best fit' line
##y = mlab.normpdf( bins, mu, sigma)
##l = plt.plot(bins, y, 'r--', linewidth=1)

#plt.xlabel('Difference in freeboard (m)')
#plt.ylabel('Occurence')
#plt.title('ALS2 - ALS1, n=125')
##plt.axis([40, 160, 0, 0.03])
#plt.grid(True)
#plt.legend()

#plt.savefig('../plots/bias')

#linear regression

mean1 = np.load('mean1.npy')
mean2 = np.load('mean2.npy')-0.05

gtri1_area = np.load('gtri1_area.npy')
gtri2_area = np.load('gtri2_area.npy')

diff = np.abs(gtri2_area-gtri1_area)/1000000
mean1 = np.ma.array(mean1,mask=diff>.005)
print mean1

als_vb1_thick = np.ma.array(als_vb1,mask=als_vb1<0.1).compressed()
als_vb2_thick = np.ma.array(als_vb2,mask=als_vb1<0.1).compressed()


#fit = np.polyfit(als_vb1_thick,als_vb2_thick,1)
fit = np.polyfit(mean1,mean2,1)
fit_fn = np.poly1d(fit) 
print fit
# fit_fn is now a function which takes in x and returns an estimate for y

#plt.plot(als_vb1_thick,als_vb2_thick, 'yo', als_vb1_thick, fit_fn(als_vb1_thick), '--k')
plt.plot(mean1,mean2, 'yo', mean1, fit_fn(mean1), '--k')

plt.savefig('../plots/bias_corr')

#





