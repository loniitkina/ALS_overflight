#! /usr/bin/python
# -*- coding: utf-8 -*-

#fix als2 specular reflections by converting nan in open water polygon to 0

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
from matplotlib.path import Path
import os

path = '../data/'
outpath = '../plots/'

#als data
#alsfile1 = '../data/ALS_20150419_EPSG3575.tiff'
alsfile2 = '../data/ALS_20150424_EPSG3575.tiff'
#open water polygon
ow2 = '../open_water_24_04'

#2nd overflight (24/04/2015)
outals2 = '../ALS_20150424'
outname2 = 'map_24_04_fix'
title2 = 'ICE-ARC overflight 24/04/2015'
llcrnrlon2=15.1; llcrnrlat2=82.72
urcrnrlon2=15.51; urcrnrlat2=82.761

#######################################333
#get the ALS data
mv =-999

#second overflight
ds = gdal.Open(alsfile2)
als2 = ds.ReadAsArray()
#this will change all nan to mv and mask them at the same time
als2 = np.ma.fix_invalid(als2, fill_value=mv)
proj = ds.GetProjection()

#this will fill all small ears of missing values, big ones will remain.
neighbors=((0,1),(0,-1),(1,0),(-1,0),(1,1),(-1,1),(1,-1),(-1,-1),
           (0,2),(0,-2),(2,0),(-2,0))

test = als2
a_copy=test.copy()
for hor_shift,vert_shift in neighbors:
    if not np.any(test.mask): break
    a_shifted=np.roll(a_copy,shift=hor_shift,axis=1)
    a_shifted=np.roll(a_shifted,shift=vert_shift,axis=0)
    idx=~a_shifted.mask*test.mask
    test[idx]=a_shifted[idx]

als2=test

#but we do not want them masked yet.
als2.mask = np.ma.nomask

#and coordinates
xx2 = np.load('xx2.npy')
yy2 = np.load('yy2.npy')

print 'got second als'

fig2    = plt.figure(figsize=(10,10))
bx      = fig2.add_subplot(111)
bx.set_title(title2, fontsize=28)

m2 = Basemap(resolution='h', epsg = 3575,
	    #ellps = 'WGS84', projection='laea', lat_0=90., lon_0=10.,
	    llcrnrlon=llcrnrlon2, llcrnrlat=llcrnrlat2,
	    urcrnrlon=urcrnrlon2, urcrnrlat=urcrnrlat2)

m2.drawcoastlines()
m2.drawparallels(np.arange(79.,90.,.01),labels=[1,0,0,0])
m2.drawmeridians(np.arange(0.,360.,.1),latmax=90.,labels=[0,0,0,1,])

m2 = Basemap(resolution='h', epsg = 3575,
	    #ellps = 'WGS84', projection='laea', lat_0=90., lon_0=10.,
	    llcrnrlon=llcrnrlon2, llcrnrlat=llcrnrlat2,
	    urcrnrlon=urcrnrlon2, urcrnrlat=urcrnrlat2)


#get the open water polygons
m2.readshapefile(ow2,'ow')

gx2, gy2 = xx2.flatten(), yy2.flatten()
points2 = np.vstack((gx2,gy2)).T

mask = np.ones_like(als2.T, dtype=bool)
mask = False
for xy in m2.ow:
    poly = shpol(xy)
    patch = PolygonPatch(poly, edgecolor='orchid', alpha=1, fill=False)
    bx.add_patch(patch)  

    #extract a mask
    px,py = poly.exterior.coords.xy
    pverts = np.vstack((px,py)).T
    path = Path(pverts)
    grid = path.contains_points(points2)
    grid = grid.reshape((xx2.shape[0],xx2.shape[1]))
    mask = np.logical_or(mask,grid)

print mask

#fix the mv in inside the survied area (specular reflections)    
als2_fixed = np.where((als2.T==mv)&(mask==True),0,als2.T)

#plot the data
im2 = m2.pcolormesh(xx2, yy2, als2_fixed, cmap=plt.cm.jet, vmin=-0.2, vmax=1.2)
plt.colorbar(im2)

#store the data
als2_fixed.dump(outals2)

fig2.tight_layout()
fig2.savefig(outpath+outname2)

