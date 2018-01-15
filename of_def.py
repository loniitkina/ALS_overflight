#! /usr/bin/python
# -*- coding: utf-8 -*-

#export virtual buoys vector layer created in QGIS to CSV files with geographical projection (WGS84, EPSG:4326)
#plot the ALS data and points

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
import osr, gdal, ogr
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon as shpol
from shapely.ops import cascaded_union
from descartes import PolygonPatch
import os
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
from datetime import date, timedelta, datetime
from of_func import *

path = '../data/'
outpath = '../plots/'

#region
llcrnrlon1=13.8; llcrnrlat1=83.1
urcrnrlon1=14.25; urcrnrlat1=83.14

llcrnrlon2=15.12; llcrnrlat2=82.725
urcrnrlon2=15.54; urcrnrlat2=82.767


#virtual buoys
vbuoys1 = '../vbuoys_1904_corr_added.csv'
vbuoys2 = '../vbuoys_2404_corr_added.csv'

vbuoys1 = '../data/vbuoys_proj_1904.csv'
vbuoys2 = '../data/vbuoys_proj_2404.csv'

#vbuoys1 = '../vbuoys_1904_new.csv'
#vbuoys2 = '../vbuoys_2404_new.csv'


#difference plots
title3 = 'ICE-ARC overflight 19-24/04/2015'
outname1 = 'def_div'
outname2 = 'def_shr'
outname3 = 'def_tdr'
########################################33

#get the virtual buoy data
data = np.loadtxt(vbuoys1, dtype=np.str,delimiter=',', skiprows=1)
latitude1 = np.array(data[:,1],dtype=float)
longitude1 = np.array(data[:,0],dtype=float)
time1 = np.array(data[:,5],dtype=float)
xo1 = np.array(data[:,2],dtype=float)
yo1 = np.array(data[:,3],dtype=float)

data = np.loadtxt(vbuoys2, dtype=np.str,delimiter=',', skiprows=1)
latitude2 = np.array(data[:,1],dtype=float)
longitude2 = np.array(data[:,0],dtype=float)
time2 = np.array(data[:,5],dtype=float)
xo2 = np.array(data[:,2],dtype=float)
yo2 = np.array(data[:,3],dtype=float)

from pyproj import Proj,Geod
g = Geod(ellps='WGS84')
epsg3575 = Proj("+init=EPSG:3575")

#get reprojected coordinates
x1v,y1v = epsg3575(longitude1,latitude1)
x2v,y2v = epsg3575(longitude2,latitude2)

#print x1[:5]
#print x2[:5]
#print x1[:5]-x2[:5]


m1 = Basemap(resolution='h', epsg = 3575,
	    #ellps = 'WGS84', projection='laea', lat_0=90., lon_0=10.,
	    llcrnrlon=llcrnrlon1, llcrnrlat=llcrnrlat1,
	    urcrnrlon=urcrnrlon1, urcrnrlat=urcrnrlat1)



m2 = Basemap(resolution='h', epsg = 3575,
	    #ellps = 'WGS84', projection='laea', lat_0=90., lon_0=10.,
	    llcrnrlon=llcrnrlon2, llcrnrlat=llcrnrlat2,
	    urcrnrlon=urcrnrlon2, urcrnrlat=urcrnrlat2)



#plot the virtual buoys
x1, y1 = m1(longitude1, latitude1)
x2, y2 = m2(longitude2, latitude2)

#x2v, y2v = m1(longitude2, latitude2)

#print x1[:5]
#print x2v[:5]
#print x1[:5]-x2v[:5]
#exit()

##get velocities from reprojected coordinates
#from pyproj import Proj,Geod
#g = Geod(ellps='WGS84')
#epsg3575 = Proj("+init=EPSG:3575")
#xp1,yp1 = epsg3575(longitude1, latitude1)
#xp2,yp2 = epsg3575(longitude2, latitude2)
#uvel = (xp1-xp2)/(5*86400)
#vvel = (yp1-yp2)/(5*86400)
#print uvel

##get velocities from original positions and times
##get time difference for each buoy
#uvel = []
#vvel = []
#for i in range(0,len(time1)):
  ##time difference
  #t1 = datetime(2015,4,19,0) + timedelta(hours=time1[i])
  #t2 = datetime(2015,4,24,0) + timedelta(hours=time2[i])
  #dt = (t2-t1).days*3600*24+(t2-t1).seconds
  ##constant time
  #dt = 5*3600*24
  ##distance
  #xd = xo1[i] - xo2[i]
  #yd = yo1[i] - yo2[i]
  ##velocity
  #uvel.append(xd/dt)
  #vvel.append(yd/dt)

#print uvel
#exit()
##print vvel
#exit()

#uvel = (x1v-x2v)/(5*86400)
#vvel = (y1v-y2v)/(5*86400)

#print xo1
#print xo2
#print xo2-xo1
#print yo2-yo1
#exit()


uvel = (xo2-xo1)/(5*86400)
vvel = (yo2-yo1)/(5*86400)


uvel = np.array(uvel)
vvel = np.array(vvel)
speed =  np.sqrt(np.mean(uvel/1000*86400)**2+np.mean(vvel/1000*86400)**2)
print 'mean speed: '+str(speed)+' km/day'

#uvel = np.ones_like(uvel)*np.mean(uvel)
#vvel = np.ones_like(vvel)*np.mean(vvel)
#print uvel
#exit()
#print vvel

#triangulate betwen the points
pts1 = np.zeros((len(x1),2))
pts2 = np.zeros((len(x1),2))
pts1[:,0]=x1; pts1[:,1]=y1
pts2[:,0]=x2; pts2[:,1]=y2

#to preserve the triangle vertices between the overflight, triangulate just the points from the 1st overflight
#then use same vertices to construct the triangles (with different coordinates) again
tri = Delaunay(pts1)
tripts1 = pts1[tri.simplices]
tripts2 = pts2[tri.simplices]

trin = len(tripts2)
print trin

dux = np.zeros(trin)
duy = np.zeros(trin)
dvx = np.zeros(trin)
dvy = np.zeros(trin)

##for each triangle
for t in range(0,trin):
  vert = tripts2[t]
  uvert = uvel[tri.simplices][t]
  vvert = vvel[tri.simplices][t]
  
  #sorting the vertices so that they are always counter-clockwise
  from scipy.spatial import ConvexHull
  hull = ConvexHull(vert)
  vert = vert[hull.vertices]
  uvert = uvert[hull.vertices]
  vvert = vvert[hull.vertices]

  dux[t],duy[t],dvx[t],dvy[t]=deformation(vert,uvert,vvert)

div = dux + dvy
div = div*1e6
shr = .5*np.sqrt((dux-dvy)**2+(duy+dvx)**2)
shr = shr*1e6
dr = np.sqrt(div**2+shr**2)

#deformation plots: divergence
fig1    = plt.figure(figsize=(10,10))
ax      = fig1.add_subplot(111)
ax.set_title(title3, fontsize=28)

m3 = Basemap(resolution='h', epsg = 3575,
	    #ellps = 'WGS84', projection='laea', lat_0=90., lon_0=10.,
	    llcrnrlon=llcrnrlon2, llcrnrlat=llcrnrlat2,
	    urcrnrlon=urcrnrlon2, urcrnrlat=urcrnrlat2)

m3.drawcoastlines()
m3.drawparallels(np.arange(79.,90.,.01),labels=[1,0,0,0])
m3.drawmeridians(np.arange(0.,360.,.1),latmax=90.,labels=[0,0,0,1,])

#virtual buoys
ax.plot(x2,y2,'o',linewidth=2,color='purple')

#triangles
patches = []
for i in range(len(tripts2)):
    patch = Polygon(tripts2[i], edgecolor='orchid', alpha=1, fill=False)
    patches.append(patch)
    
p = PatchCollection(patches, cmap=plt.cm.bwr, alpha=0.4)
p.set_array(np.array(div))
ax.add_collection(p)
p.set_clim([-3, 3])
cbar = plt.colorbar(p)
cbar.set_label(r'Divergence (10$^6$s$^{-1}$)',size=16)

fig1.tight_layout()
fig1.savefig(outpath+outname1)

#deformation plots: shear
fig2    = plt.figure(figsize=(10,10))
bx      = fig2.add_subplot(111)
bx.set_title(title3, fontsize=28)

m3 = Basemap(resolution='h', epsg = 3575,
	    #ellps = 'WGS84', projection='laea', lat_0=90., lon_0=10.,
	    llcrnrlon=llcrnrlon2, llcrnrlat=llcrnrlat2,
	    urcrnrlon=urcrnrlon2, urcrnrlat=urcrnrlat2)

m3.drawcoastlines()
m3.drawparallels(np.arange(79.,90.,.01),labels=[1,0,0,0])
m3.drawmeridians(np.arange(0.,360.,.1),latmax=90.,labels=[0,0,0,1,])

#virtual buoys
bx.plot(x2,y2,'o',linewidth=2,color='purple')

#triangles
patches = []
for i in range(len(tripts2)):
    patch = Polygon(tripts2[i], edgecolor='orchid', alpha=1, fill=False)
    patches.append(patch)
    
p = PatchCollection(patches, cmap=plt.cm.Reds, alpha=0.4)
p.set_array(np.array(shr))
bx.add_collection(p)
p.set_clim([0, 3])
cbar = plt.colorbar(p)
cbar.set_label(r'Shear (10$^6$s$^{-1}$)',size=16)

fig2.tight_layout()
fig2.savefig(outpath+outname2)

#deformation plots: total deformation
fig3    = plt.figure(figsize=(10,10))
cx      = fig3.add_subplot(111)
cx.set_title(title3, fontsize=28)

m3 = Basemap(resolution='h', epsg = 3575,
	    #ellps = 'WGS84', projection='laea', lat_0=90., lon_0=10.,
	    llcrnrlon=llcrnrlon2, llcrnrlat=llcrnrlat2,
	    urcrnrlon=urcrnrlon2, urcrnrlat=urcrnrlat2)

m3.drawcoastlines()
m3.drawparallels(np.arange(79.,90.,.01),labels=[1,0,0,0])
m3.drawmeridians(np.arange(0.,360.,.1),latmax=90.,labels=[0,0,0,1,])

#virtual buoys
cx.plot(x2,y2,'o',linewidth=2,color='purple')

#triangles
patches = []
for i in range(len(tripts2)):
    patch = Polygon(tripts2[i], edgecolor='orchid', alpha=1, fill=False)
    patches.append(patch)
    
p = PatchCollection(patches, cmap=plt.cm.Reds, alpha=0.4)
p.set_array(np.array(dr))
cx.add_collection(p)
p.set_clim([0, 3])
cbar = plt.colorbar(p)
cbar.set_label(r'Total deformation (10$^6$s$^{-1}$)',size=16)

fig3.tight_layout()
fig3.savefig(outpath+outname3)

#np.save('div',div)
#np.save('shr',shr)