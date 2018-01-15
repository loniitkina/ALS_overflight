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
from matplotlib.patches import Polygon
from shapely.ops import cascaded_union
import csv
from mpl_toolkits.axes_grid.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import date, timedelta, datetime, time
import osr, gdal, ogr
#from osgeo import gdal, osr
from of_func import *

path = '../data/frm_DTU_020816/ALS/'
#path = '../data/frm_DTU_040417/'
outpath = '../plots/'
outpath_data = '../data/'

#ALS1
alsfile = path+'ALS_20150419_1032_5x5_corr_sh_dc.frs'
#alsfile = path+'ALS_20150419_1032_5x5_raw_dc.frs'
#open water polygon
ow = '../open_water_19_04'
vbuoys = '../vbuoys_1904_new.csv'
vbuoys = '../vbuoys_1904_thinned.csv'
vbuoys = '../vbuoys_1904_revision.csv'
gtri = np.load('gtri1.npy')
vbuoys_out = 'vbuoys_proj_1904.csv'
reg = 'als1_zoom'
name_als = 'als19_fb_zoom.png'
name_time = 'als19_tm_zoom.png'
slon = 14.16;slat = 83.097
lance_lon = 14.064979223402144;lance_lat = 83.122550903451426
name_tiff = 'als19_fb_zoom.tiff'
name_dump = 'als19_fb'
name_dumpx = 'als19_x'
name_dumpy = 'als19_y'
name_dumplon = 'als19_lon'
name_dumplat = 'als19_lat'
title = 'a'

#ALS2
alsfile = path+'ALS_20150424_0955_5x5_corr_sh_dc.frs'
#alsfile = path+'ALS_20150424_0955_5x5_raw_dc.frs'
#open water polygon
ow = '../open_water_24_04'
waves = '../waves_24_04'
vbuoys = '../vbuoys_2404_new.csv'
vbuoys = '../vbuoys_2404_thinned.csv'
vbuoys = '../vbuoys_2404_revision.csv'
gtri = np.load('gtri2.npy')
vbuoys_out = 'vbuoys_proj_2404.csv'
reg = 'als2_zoom'
name_als = 'als24_fb_zoom.png'
name_time = 'als24_tm_zoom.png'
slon = 15.44;slat = 82.718
#Lance
lance_lon = 15.331769876266121;lance_lat = 82.743808335535192
name_tiff = 'als24_fb_zoom.tiff'
name_dump = 'als24_fb'
name_dumpx = 'als24_x'
name_dumpy = 'als24_y'
name_dumplon = 'als24_lon'
name_dumplat = 'als24_lat'
title = 'b'

#pick a reasonable min value to get the majority of freeboars above zero in both overflights
minval = -.15


#read als data
als = np.array(getColumn(alsfile,11),dtype=float)
#als_elev = np.array(getColumn(alsfile,12),dtype=float)
#als=als_elev	#quick check if raw elevations are better than processed FBs!
tm = np.array(getColumn(alsfile,0),dtype=np.datetime64).astype(datetime)
tm = tm.tolist()
tm_hr = np.array(map(lambda x: x.hour + (x.minute)/60. + (x.second)/3600., tm))
lons = np.array(getColumn(alsfile,4),dtype=float)
lats = np.array(getColumn(alsfile,2),dtype=float)
xshift = np.array(getColumn(alsfile,10),dtype=float)
yshift = np.array(getColumn(alsfile,9),dtype=float)

#get the virtual buoy data
data = np.loadtxt(vbuoys, dtype=np.str,delimiter=',', skiprows=1)
latitude1 = np.ma.array(data[:,1],mask=data[:,1]=='0',dtype=float)
longitude1 = np.ma.array(data[:,0],mask=data[:,2]=='0',dtype=float)

#define projection and interpolate to regular grid
area_def = pr.utils.load_area('area.cfg', reg)
swath_def = geometry.SwathDefinition(lons=lons, lats=lats)
als_map = kd_tree.resample_gauss(swath_def, als,area_def, radius_of_influence=10, sigmas=5)
#als_map = kd_tree.resample_nearest(swath_def, als,area_def, radius_of_influence=10)
#time_map = kd_tree.resample_gauss(swath_def, time,area_def, radius_of_influence=10, sigmas=5)
time_map = kd_tree.resample_nearest(swath_def, tm_hr,area_def, radius_of_influence=10)
xshift_map = kd_tree.resample_nearest(swath_def, xshift,area_def, radius_of_influence=10)
yshift_map = kd_tree.resample_nearest(swath_def, yshift,area_def, radius_of_influence=10)

#mask irrelevant and extreme data
#mask_bad = (time_map<10)|(time_map>14)|(als_map<-1)|(als_map>2)
#mask_bad = (als_map==0)
#als_map = np.ma.array(als_map, mask=mask_bad)
#time_map = np.ma.array(time_map, mask=mask_bad)


##plotting
#pr.plot.save_quicklook(outpath+'als14_quick.png', area_def, result, label='ALS freeboard (m)')
#ALS data
fig1    = plt.figure(figsize=(10,10))
cx      = fig1.add_subplot(111)
cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

bmap1 = pr.plot.area_def2basemap(area_def)
print als_map.shape
x,y = np.meshgrid(np.arange(bmap1.xmin,bmap1.xmax-1,5),np.arange(bmap1.ymax,bmap1.ymin+1,-5))
print x.shape
if title =='b':
  #fix specular reflections
  #get the open water polygons
  bmap1.readshapefile(ow,'ow', drawbounds=False)

  x,y = np.meshgrid(np.arange(bmap1.xmin,bmap1.xmax-1,5),np.arange(bmap1.ymax,bmap1.ymin+1,-5))
  gx2, gy2 = x.flatten(), y.flatten()
  points2 = np.vstack((gx2,gy2)).T

  mask = np.ones_like(als_map, dtype=bool)
  mask = False
  for xy in bmap1.ow:
      poly = shpol(xy)

      #extract a mask
      px,py = poly.exterior.coords.xy
      pverts = np.vstack((px,py)).T
      path = Path(pverts)
      grid = path.contains_points(points2)
      grid = grid.reshape((x.shape[0],x.shape[1]))
      mask = np.logical_or(mask,grid)

  #print mask
  #mask[0,0] = True
  #print np.min(mask), np.max(mask)
  #print mask.shape, grid.shape, als_map.shape

  #fix the mv in inside the survied area (specular reflections) and add the 'lead bias'
  als_map = np.where((als_map==0)&(mask==True),minval,als_map)
  
  #fix the waves in open water
  bmap1.readshapefile(waves,'waves', drawbounds=False)

  x,y = np.meshgrid(np.arange(bmap1.xmin,bmap1.xmax-1,5),np.arange(bmap1.ymax,bmap1.ymin+1,-5))
  gx2, gy2 = x.flatten(), y.flatten()
  points2 = np.vstack((gx2,gy2)).T

  mask = np.ones_like(als_map, dtype=bool)
  mask = False
  for xy in bmap1.waves:
      poly = shpol(xy)

      #extract a mask
      px,py = poly.exterior.coords.xy
      pverts = np.vstack((px,py)).T
      path = Path(pverts)
      grid = path.contains_points(points2)
      grid = grid.reshape((x.shape[0],x.shape[1]))
      mask = np.logical_or(mask,grid)
  als_map = np.where(mask==True,minval,als_map)

#estimate bias
als_map = np.where((als_map<minval),minval,als_map)
lead = minval*-1
  
mv=-999  
als_map = np.where((als_map==0)&(time_map==0),mv,als_map+lead)
mask_bad = (als_map==mv)
als_map = np.ma.array(als_map, mask=mask_bad)
lancex, lancey = bmap1(lance_lon, lance_lat)
bmap1.plot(lancex,lancey,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

col = bmap1.imshow(als_map, origin='upper', vmin=0, vmax=1)
x1, y1 = bmap1(longitude1, latitude1)
bmap1.plot(x1,y1,'o',linewidth=2,color='purple',markersize=8)
#import pdb; pdb.set_trace()	#debug
#triangles
patches = []
for i in range(len(gtri)):
    patch = Polygon(gtri[i], edgecolor='orchid',linewidth=2, alpha=1, fill=False)
    poly = shpol(gtri[i])
    patches.append(poly)
    cx.add_patch(patch)

###unified triangles   
#group = cascaded_union(patches)
#patch = PolygonPatch(group, edgecolor='k', alpha=1, fill=False)
#cx.add_patch(patch)

## create an axes on the right side of ax. The width of cax will be 5%
## of ax and the padding between cax and ax will be fixed at 0.05 inch.
#divider = make_axes_locatable(cx)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#cbar = bmap1.colorbar(col, cax=cax)
cbar = bmap1.colorbar(location='bottom')
cbar.set_label(r'Freeboard (m)',size=16)

bmap1.drawmapscale(slon, slat, 10, 90, 1000, units='m', barstyle='fancy',fontsize=14)
#North arrow
xa, ya = bmap1(slon, slat)
plt.arrow(xa+800,ya,0,75,fc="k", ec="k", linewidth = 4, head_width=200, head_length=200,overhang=.5)
plt.text(xa+760,ya-150,'N')
fig1.savefig(outpath+name_als, bbox_inches='tight')

#time
fig2    = plt.figure(figsize=(10,10))
cx      = fig2.add_subplot(111)
cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

bmap2 = pr.plot.area_def2basemap(area_def)
col = bmap2.imshow(time_map, origin='upper')
x1, y1 = bmap2(longitude1, latitude1)
bmap2.plot(x1,y1,'o',linewidth=2,color='purple')
bmap2.colorbar()
fig2.savefig(outpath+name_time, bbox_inches='tight')


#get coordinates, als and time for all the buoys
#get projected coordinates for the whole region
x,y = area_def.get_proj_coords()
lons,lats = area_def.get_lonlats()
#this is how i get indexes of buoy coordinates in the array of projected coordinates (x,y)
biy,bix = area_def.get_xy_from_lonlat(longitude1,latitude1)
#this gives me projected x coordinate of the first buoy
#print x[bix[0],biy[0]]
#find nearest als value and time value for all the buoys (should have same index!)
bx = []
by = []
bals = []
btm = []
bxsh = []
bysh = []
for i in range(0,len(bix)):
  #print bix[i],biy[i]
  bx.append(x[bix[i],biy[i]])
  by.append(y[bix[i],biy[i]])
  bals.append(als_map[bix[i],biy[i]])
  btm.append(time_map[bix[i],biy[i]])
  bxsh.append(xshift_map[bix[i],biy[i]])  
  bysh.append(yshift_map[bix[i],biy[i]])

##calculate original position
bxo = np.array(bx) - np.array(bxsh)
byo = np.array(by) - np.array(bysh)

#save a text file (table)
tt = [longitude1,latitude1,bx,by,bals,btm,bxo,byo]
table = zip(*tt)
#print table[0]

outname = outpath_data+vbuoys_out
with open(outname, 'wb') as f:
  #header
  f.write(b'lon,lat,x(m),y(m),ALS freeboard (m),time (decimal hour)\n')
  np.savetxt(f, table, fmt="%s", delimiter=",")


#save all arrays
als_map.mask = np.ma.nomask
als_map.dump(outpath_data+name_dump)
x.dump(outpath_data+name_dumpx)
y.dump(outpath_data+name_dumpy)
lons.dump(outpath_data+name_dumplon)
lats.dump(outpath_data+name_dumplat)


exit()

def save_geotiff(raster_array, area_def, export_path):
	
	# set manually the number of bands to one, because we know you have only one layer
	bands_number = 1
	
	# Create gdal geotiff driver
        gtiff_driver = gdal.GetDriverByName('GTiff')

	# Pick up the numbers format, use gdal.GDT_Float32 for floats
        gtiff_format = gdal.GDT_Float64
        gtiff_options=["COMPRESS=LZW", "PREDICTOR=2", "TILED=YES"]
        gtiff_options = ["COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6", "INTERLEAVE=BAND"]
        gtiff_options = []
        
	# Create output file (empty at this point) and define projection and extent parameters
	gtiff_dataset = gtiff_driver.Create(export_path,
                                             int(area_def.x_size),
                                             int(area_def.y_size),
                                             bands_number,
                                             gtiff_format,
                                             gtiff_options)

	# Define area extent for the Geotiff dataset        
	geometry_list = (area_def.area_extent[0],
                         area_def.pixel_size_x,
                         0,
                         area_def.area_extent[3],
                         0,
                         area_def.pixel_size_y * -1)

	# Set projection parameters
        gtiff_dataset.SetGeoTransform(geometry_list)
        srs = osr.SpatialReference()
        srs.ImportFromProj4(area_def.proj4_string.encode('ascii'))
        gtiff_dataset.SetProjection(srs.ExportToWkt())
	
	# Get the empty band from the dataset, so gdal knows where to write the data	
	gtiff_band = gtiff_dataset.GetRasterBand(1)
	
	# Write the layer (your raster array) data into the geotiff dataset
	gtiff_band.WriteArray(raster_array)

	gtiff_dataset = None

geotiff_file = outpath_data+name_tiff
save_geotiff(als_map, area_def, geotiff_file)
