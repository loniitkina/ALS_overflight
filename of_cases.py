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
import pyresample as pr
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon as shpol
from shapely.ops import cascaded_union
from descartes import PolygonPatch
import os
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.path import Path
from scipy.stats import mode
from of_func import *

path = '../data/'
outpath = '../plots/'
outpath_tri = '../plots/pdf_tri/'

mv = -999.

#als data
als1 = np.load(path+'als19_fb')
lons1 = np.load(path+'als19_lon')
lats1 = np.load(path+'als19_lat')

als2 = np.load(path+'als24_fb')
lons2 = np.load(path+'als24_lon')
lats2 = np.load(path+'als24_lat')

#virtual buoys
vbuoys1 = '../vbuoys_1904_thinned.csv'
vbuoys2 = '../vbuoys_2404_thinned.csv'
#######################################33
outname_legend = 'case_legend'
date1 = '19.04.2015'
date2 = '24.04.2015'

##rafting
#outname1 = 'case_raft'
#title1 = 'a'
#reg1 = 'als1_raft'
#reg2 = 'als2_raft'

##ridging
#outname1 = 'case_ridge'
#title1 = 'b'
#reg1 = 'als1_ridge'
#reg2 = 'als2_ridge'

#new lead
outname1 = 'case_lead'
title1 = 'c'
reg1 = 'als1_lead'
reg2 = 'als2_lead'
##########################################33

#get the virtual buoy data
data = np.loadtxt(vbuoys1, dtype=np.str,delimiter=',', skiprows=1)
latitude1 = np.ma.array(data[:,1],dtype=float)
longitude1 = np.ma.array(data[:,0],dtype=float)

data = np.loadtxt(vbuoys2, dtype=np.str,delimiter=',', skiprows=1)
latitude2 = np.ma.array(data[:,1],dtype=float)
longitude2 = np.ma.array(data[:,0],dtype=float)

#plotting part
area_def1 = pr.utils.load_area('area.cfg', reg1)
area_def2 = pr.utils.load_area('area.cfg', reg2)
m2 = pr.plot.area_def2basemap(area_def2)

#first plot
fig1   = plt.figure(figsize=(10,10))
plt.subplots_adjust(hspace = .001)
ax     = fig1.add_subplot(211)
ax.plot(.05, .9, 'w.', markersize=110, transform=ax.transAxes, markeredgecolor='k', markeredgewidth=2)
ax.text(.05, .9, title1, ha='center', va='center', transform=ax.transAxes, fontdict={'color':'k','size':55})
props = dict(boxstyle='round', facecolor='w', alpha=1)
ax.text(.87, .93, date1, ha='center', va='center', transform=ax.transAxes, fontdict={'color':'k','size':24},bbox=props)
m1 = pr.plot.area_def2basemap(area_def1)

#als data
xx1,yy1 = m1(lons1,lats1)
xx2,yy2 = m2(lons2,lats2)
als1 = np.ma.array(als1,mask=als1==-999)
im1 = m1.pcolormesh(xx1, yy1, als1, cmap=plt.cm.jet, vmin=0, vmax=1)
als1.mask = np.ma.nomask

# Then, I wanted to have a scale on the plot :
asb =  AnchoredSizeBar(ax.transData,100.,"100 m", loc=4,pad=.5, borderpad=.25, sep=5,frameon=True)
ax.add_artist(asb) # and finaly, we add the scale to the inset axis

#second plot
bx     = fig1.add_subplot(212)
bx.text(.87, .93, date2, ha='center', va='center', transform=bx.transAxes, fontdict={'color':'k','size':24},bbox=props)
m2 = pr.plot.area_def2basemap(area_def2)

als2 = np.ma.array(als2,mask=als2==-999)
im2 = m2.pcolormesh(xx2, yy2, als2, cmap=plt.cm.jet, vmin=0, vmax=1)
als2.mask = np.ma.nomask

# Then, I wanted to have a scale on the plot :
asb =  AnchoredSizeBar(bx.transData,100.,"100 m", loc=4,pad=.5, borderpad=.25, sep=5,frameon=True)
bx.add_artist(asb) # and finaly, we add the scale to the inset axis

#virtual buoys
x1, y1 = m1(longitude1, latitude1)
x2, y2 = m2(longitude2, latitude2)
ax.plot(x1,y1,'o',markersize=10,color='purple')
bx.plot(x2,y2,'o',markersize=10,color='purple')


#triangulation
pts1 = np.zeros((len(x1),2))
pts2 = np.zeros((len(x1),2))
pts1[:,0]=x1; pts1[:,1]=y1
pts2[:,0]=x2; pts2[:,1]=y2
#to preserve the triangle vertices between the overflight, triangulate just the points form the 1st overflight
#then use same vertices to construct the triangles (with different coordinates) again
tri = Delaunay(pts1)
tripts1 = pts1[tri.simplices]
tripts2 = pts2[tri.simplices]

gx1, gy1 = xx1.flatten(), yy1.flatten()
points1 = np.vstack((gx1,gy1)).T
gx2, gy2 = xx2.flatten(), yy2.flatten()
points2 = np.vstack((gx2,gy2)).T

#triangles
patches = []
for i in range(len(tripts1)):
    #check if there are any missing values in the triangle and append only triangles that cover fully data-covered areas
    path1 = Path(tripts1[i])
    grid1 = path1.contains_points(points1)
    grid1 = grid1.reshape((xx1.shape[0],xx1.shape[1]))
    
    path2 = Path(tripts2[i])
    grid2 = path2.contains_points(points2)
    grid2 = grid2.reshape((xx2.shape[0],xx2.shape[1]))
    
    #use inverted grid as mask for the als array
    tmp1 = np.ma.array(als1,mask=~grid1)
    tmp2 = np.ma.array(als2,mask=~grid2)

    if (mv in np.ma.compressed(tmp1)) or (mv in np.ma.compressed(tmp2)):
      continue
    else:
      patch = Polygon(tripts1[i], edgecolor='orchid',linewidth=4, alpha=.8, fill=False)
      ax.add_patch(patch)
      patch = Polygon(tripts2[i], edgecolor='orchid',linewidth=4, alpha=.8, fill=False)
      bx.add_patch(patch)

fig1.tight_layout()
fig1.savefig(outpath+outname1)


# Just colorbar
fig2 = plt.figure(figsize=(6,1.2))
#ax1 = fig2.add_subplot(121)
ax1 = fig2.add_axes([.01, .5, .95, .4])

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
cmap = plt.cm.jet
import matplotlib as mpl
norm = mpl.colors.Normalize(vmin=0, vmax=1)

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('Freeboard (m)')

#fig2.tight_layout()
fig2.savefig(outpath+outname_legend)
