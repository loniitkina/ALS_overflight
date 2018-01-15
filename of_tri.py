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

#als data
als1 = np.load(path+'als19_fb')
#xx1 = np.load(path+'als19_x')
#yy1 = np.load(path+'als19_y')
lons1 = np.load(path+'als19_lon')
lats1 = np.load(path+'als19_lat')
reg1 = 'als1_zoom'

als2 = np.load(path+'als24_fb')
#xx2 = np.load(path+'als24_x')
#yy2 = np.load(path+'als24_y')
lons2 = np.load(path+'als24_lon')
lats2 = np.load(path+'als24_lat')
reg2 = 'als2_zoom'

#virtual buoys
vbuoys1 = '../vbuoys_1904_new.csv'
vbuoys2 = '../vbuoys_2404_new.csv'

vbuoys1 = '../vbuoys_1904_thinned.csv'
vbuoys2 = '../vbuoys_2404_thinned.csv'
#######################################33
#1st overflight (19/04/2015)
outals1 = '../ALS_20150419'
outname1 = 'map_19_04'
title1 = 'a'

#2nd overflight (24/04/2015)
outals2 = '../ALS_20150424'
outname2 = 'map_24_04'
title2 = 'b'
########################################33

#get the virtual buoy data
data = np.loadtxt(vbuoys1, dtype=np.str,delimiter=',', skiprows=1)
latitude1 = np.ma.array(data[:,1],dtype=float)
longitude1 = np.ma.array(data[:,0],dtype=float)

data = np.loadtxt(vbuoys2, dtype=np.str,delimiter=',', skiprows=1)
latitude2 = np.ma.array(data[:,1],dtype=float)
longitude2 = np.ma.array(data[:,0],dtype=float)

#get the ALS data
mv =-999

#interpolate some holes
als1 = np.ma.array(als1,mask=als1==mv)
neighbors=((0,1),(0,-1),(1,0),(-1,0),(1,1),(-1,1),(1,-1),(-1,-1),
           (0,2),(0,-2),(2,0),(-2,0),(-2,1),(-2,-2),(-2,-1),(-2,1),(-2,2),(-1,2),(2,1),(2,2),(1,2),(2,-1),(2,-2),(1,-2))

a_copy=als1.copy()
for hor_shift,vert_shift in neighbors:
  if not np.any(als1.mask): break
  a_shifted=np.roll(a_copy,shift=hor_shift,axis=1)
  a_shifted=np.roll(a_shifted,shift=vert_shift,axis=0)
  idx=~a_shifted.mask*als1.mask
  als1[idx]=a_shifted[idx]

als1.mask = np.ma.nomask 

als2 = np.ma.array(als2,mask=als2==mv)
a_copy=als2.copy()
for hor_shift,vert_shift in neighbors:
    if not np.any(als2.mask): break
    a_shifted=np.roll(a_copy,shift=hor_shift,axis=1)
    a_shifted=np.roll(a_shifted,shift=vert_shift,axis=0)
    idx=~a_shifted.mask*als2.mask
    als2[idx]=a_shifted[idx]

als2.mask = np.ma.nomask 

area_def1 = pr.utils.load_area('area.cfg', reg1)
m1 = pr.plot.area_def2basemap(area_def1)
area_def2 = pr.utils.load_area('area.cfg', reg2)
m2 = pr.plot.area_def2basemap(area_def2)

xx1,yy1 = m1(lons1,lats1)
xx2,yy2 = m2(lons2,lats2)

#plot the virtual buoys
x1, y1 = m1(longitude1, latitude1)
x2, y2 = m2(longitude2, latitude2)
x2v, y2v = m1(longitude2, latitude2)

print len(x1)
print len(x2)
	
#triangulate betwen the points
pts1 = np.zeros((len(x1),2))
pts2 = np.zeros((len(x1),2))
pts1[:,0]=x1; pts1[:,1]=y1
pts2[:,0]=x2; pts2[:,1]=y2

#to preserve the triangle vertices between the overflight, triangulate just the points form the 1st overflight
#then use same vertices to construct the triangles (with different coordinates) again
tri = Delaunay(pts1)
tripts1 = pts1[tri.simplices]
tripts2 = pts2[tri.simplices]

uvel = (x2v-x1)/(5*86400)
vvel = (y2v-y1)/(5*86400)
upts = uvel[tri.simplices]
vpts = vvel[tri.simplices]
print 'done triangulation'

#make a plot of buoys with the velocities (color of the marker depending on the velocity)
fig1    = plt.figure(figsize=(10,10))
cx      = fig1.add_subplot(111)
cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
cx.text(.05, .95, 'f', ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

m1.drawcoastlines()

#virtual buoys
speed = np.sqrt(uvel**2+vvel**2)*100	#m/s to cm/s
sc = cx.scatter(x2,y2,s=500,c=speed,cmap=plt.cm.Reds)
# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(cx)
cax = divider.append_axes("bottom", size="5%", pad=0.1)
cbar = plt.colorbar(sc, cax=cax, ticks=np.arange(10.44,10.68,.02),orientation='horizontal')
cbar.set_label('speed (cm/s)',size=16)

fig1.savefig(outpath+'map_vel_speed',bbox_inches='tight')

#filter out only good triangles
gtri1=[]
gtri2=[]
gtri1_nc=[]
gtri2_nc=[]
gtri1_bc=[]
gtri2_bc=[]
gtri1_bc_pos=[]
gtri2_bc_pos=[]
gtri1_bc_neg=[]
gtri2_bc_neg=[]

gtri1_area=[]
gtri2_area=[]
gtri1_nc_area=[]
gtri2_nc_area=[]
gtri1_bc_area=[]
gtri2_bc_area=[]
gtri1_bc_pos_area=[]
gtri2_bc_pos_area=[]
gtri1_bc_neg_area=[]
gtri2_bc_neg_area=[]

mean1=[]
mean2=[]
mean1_bc_pos=[]
mean2_bc_pos=[]
mean1_bc_neg=[]
mean2_bc_neg=[]

gtri1_nc_fb=[]
gtri2_nc_fb=[]
gtri1_bc_pos_fb=[]
gtri2_bc_pos_fb=[]
gtri1_bc_neg_fb=[]
gtri2_bc_neg_fb=[]

mode1=[]
mode2=[]

bias = []

dux = []
duy = []
dvx = []
dvy = []
gtri_minang = []
gtri_nc_div=[]
gtri_bc_pos_div=[]
gtri_bc_neg_div=[]

#save just simplices/triangles ids
#alternative masking per triangle
gx1, gy1 = xx1.flatten(), yy1.flatten()
points1 = np.vstack((gx1,gy1)).T
gx2, gy2 = xx2.flatten(), yy2.flatten()
points2 = np.vstack((gx2,gy2)).T

offset = np.ones(len(tripts1))*.03
offset[1]=.07
offset[3]=.08
offset[4]=.09
offset[5]=-.05
offset[6]=-.05
offset[12]=.04
offset[13]=.05	#from .04
offset[14]=.04
offset[15]=.02
offset[17]=.02	#from .04
offset[20]=-.03
offset[21]=-.1
offset[22]=-.1
offset[27]=.01
offset[28]=.0
offset[29]=.01
offset[30]=.01
offset[31]=.0
offset[35]=.02
offset[36]=-.01
offset[40]=.02
offset[44]=.05
offset[45]=.02
offset[46]=-.09
offset[47]=.05
offset[48]=.1
offset[51]=.01
offset[52]=.0
offset[58]=-.03
offset[59]=-.03
offset[60]=-.08
offset[62]=-.07
offset[63]=-.07
offset[64]=.04
offset[65]=.04
offset[68]=.04
offset[71]=-.04
offset[73]=.0
offset[79]=.08
offset[80]=.12
offset[81]=.13
offset[85]=.02
offset[86]=-.01
offset[88]=.08
offset[89]=.08
offset[90]=.02
offset[91]=.15
offset[92]=.01
offset[93]=.1
offset[94]=.04
offset[95]=.06
offset[96]=.03
offset[97]=.05
offset[98]=.06
offset[99]=-.02
offset[101]=.01
offset[102]=.01
offset[103]=.03
offset[104]=.02
offset[105]=.03
offset[106]=.03
offset[109]=-.01
offset[110]=-.13 #from -.06
offset[111]=-.08 #from -.03
offset[112]=.01
offset[113]=-.02
offset[114]=-.1 #from -.18
offset[116]=-.15 #from -.18
offset[118]=-.06 #from -.05
offset[120]=-.05 #from -.03
offset[121]=.02
offset[123]=.04
offset[124]=.04
offset[125]=.02
offset[126]=.01
offset[127]=.04
offset[128]=.06
offset[129]=.06
offset[130]=-.05
offset[131]=.0
offset[132]=-.06
offset[133]=-.03
offset[134]=-.06
offset[135]=-.03
offset[136]=-.09 #from -.1
offset[137]=.07
offset[138]=.03
offset[139]=.11
offset[140]=.08
offset[141]=.1
offset[142]=.07
offset[143]=.13
offset[144]=.11
offset[145]=.02
offset[146]=-.02
offset[147]=.02
offset[148]=.01
offset[152]=.05
offset[153]=.02
offset[154]=.15 #from .18
offset[155]=.03
offset[156]=.1
offset[157]=.03
offset[158]=.01
offset[159]=.03
offset[160]=.08
offset[161]=.07
offset[162]=.09
offset[163]=.1
offset[164]=.09
offset[166]=.01
offset[167]=.03
offset[168]=.08
offset[169]=.11
offset[170]=.05
offset[171]=.04
offset[172]=.03
offset[173]=.02
offset[174]=.04
offset[175]=.03
offset[176]=.09
offset[177]=.0
offset[178]=.04
offset[179]=.12 #from .11
offset[180]=.0 #from .02
offset[181]=.06
offset[182]=.03
offset[183]=-.02 #from .02
offset[184]=-.13
offset[185]=.05
offset[186]=.03
offset[187]=.03
offset[188]=.07
offset[189]=.09
offset[190]=.05
offset[191]=.13
offset[192]=.05
offset[193]=.04
offset[194]=.05
offset[195]=.03
offset[196]=.02
offset[197]=.03
offset[198]=.04
offset[199]=.07
offset[200]=.08
offset[201]=.08
offset[202]=.11
offset[203]=.07
offset[204]=.04
offset[205]=.07
offset[206]=.07
offset[207]=.05
offset[210]=.03
offset[211]=.02

print len(tripts1)
for t in range(0,len(tripts1)):
    vert1 = np.asarray(tripts1[t])
    vert2 = np.asarray(tripts2[t])
    uvert = upts[t]
    vvert = vpts[t]

    #sorting the vertices so that they are always counter-clockwise
    hull = ConvexHull(vert1)
    vert1 = vert1[hull.vertices]
    
    hull = ConvexHull(vert2)
    vert2 = vert2[hull.vertices]
    uvert = uvert[hull.vertices]
    vvert = vvert[hull.vertices]

    #check if there are any missing values in the triangle and append only triangles that cover fully data-covered areas
    path1 = Path(vert1)
    grid1 = path1.contains_points(points1)
    grid1 = grid1.reshape((xx1.shape[0],xx1.shape[1]))
    
    path2 = Path(vert2)
    grid2 = path2.contains_points(points2)
    grid2 = grid2.reshape((xx2.shape[0],xx2.shape[1]))
    
    #use inverted grid as mask for the als array
    tmp1 = np.ma.array(als1,mask=~grid1)
    tmp2 = np.ma.array(als2,mask=~grid2)
    
    #print np.ma.compressed(tmp1)
     
    if (mv in np.ma.compressed(tmp1)) or (mv in np.ma.compressed(tmp2)):
      print 'ommiting triangle: ', t
    else:
      gtri1.append(vert1)
      gtri2.append(vert2)
      
      #calculate deformation
      a,b,c,d,minang=deformation(vert2,uvert,vvert)
      dux.append(a)
      duy.append(b)
      dvx.append(c)
      dvy.append(d)
      ddd = a + d
      sss = .5*np.sqrt((a-d)**2+(b+c)**2)
      
      #store minangle
      gtri_minang.append(minang)

      #calculate area
      area1 = .5* (vert1[0,0]*vert1[1,1] - vert1[0,1]*vert1[1,0] + vert1[1,0]*vert1[2,1] - vert1[1,1]*vert1[2,0] + vert1[2,0]*vert1[0,1] - vert1[2,1]*vert1[0,0])
      area2 = .5* (vert2[0,0]*vert2[1,1] - vert2[0,1]*vert2[1,0] + vert2[1,0]*vert2[2,1] - vert2[1,1]*vert2[2,0] + vert2[2,0]*vert2[0,1] - vert2[2,1]*vert2[0,0])
      gtri1_area.append(area1)
      gtri2_area.append(area2)
           
      #check is there is an offest in the thick ice freeboards (no growth possible, no deformation there)
      ttmp1 = np.ma.compressed(tmp1)#+ .03
      ttmp2 = np.ma.compressed(tmp2)
      
      from scipy import stats
      #binval=10
      #bins = np.linspace(.2, 1,num=binval)
      #count,lowerlimit,binsize,extra= stats.histogram(ttmp1,binval,defaultlimits=(.2, 1))
      #mindx1 =  np.argmax(count)
      #count,lowerlimit,binsize,extra= stats.histogram(ttmp2,binval,defaultlimits=(.2, 1))
      #mindx2 =  np.argmax(count)
      #offset = bins[mindx2] - bins[mindx1]
      #if  np.abs(ddd*1e6) < .3:
	#ttmp1 = ttmp1 + offset
      ##in case of strong divergence or convergence do something special...
      #else:
	#ttmp1 = ttmp1 + .03
      ttmp1 = ttmp1 + offset[t]
      
      bias.append(offset[t])
      
      #get mean and mode
      me1 = np.mean(ttmp1)
      me2 = np.mean(ttmp2)
      mean1.append(me1)
      mean2.append(me2)
      
      #primary mode
      binval=200
      bins = np.linspace(-.2, 2,num=binval)
      count,lowerlimit,binsize,extra= stats.histogram(ttmp1,binval,defaultlimits=(-.2, 2))
      mindx1 =  np.argmax(count)
      count,lowerlimit,binsize,extra= stats.histogram(ttmp2,binval,defaultlimits=(-.2, 2))
      mindx2 =  np.argmax(count)
      m1 = bins[mindx2]
      m2 = bins[mindx1]
      mode1.append(m1)
      mode2.append(m2)
      
      #diff = (area2-area1)/1000000
      #diff_r = diff/(area1/1000000)
      ##pdf plot
      #bins = 120
      #fig1= plt.figure(figsize=(6,6))
      #ax = fig1.add_subplot(111)
      #ax.hist(ttmp1, bins, normed=False, histtype='stepfilled', alpha=0.4, lw=0, facecolor='red',  range=[-.2, 2])
      #ax.hist(ttmp2, bins, normed=False, histtype='stepfilled', alpha=0.4, lw=0, facecolor='blue',  range=[-.2, 2])
      #ax.text(.05,.95,'mode 1: '+str(m1), transform = ax.transAxes)
      #ax.text(.05,.9,'mode 2: '+str(m2), transform = ax.transAxes)
      #ax.text(.05,.85,'mean 1: '+str(me1), transform = ax.transAxes)
      #ax.text(.05,.8,'mean 2: '+str(me2), transform = ax.transAxes)
      #ax.text(.05,.75,'div : '+str(np.round(ddd*1e6,2)), transform = ax.transAxes)
      #ax.text(.05,.7,'shr : '+str(np.round(sss*1e6,2)), transform = ax.transAxes)
      #ax.text(.05,.65,'offset: '+str(offset[t]), transform = ax.transAxes)
      #ax.text(.05,.6,'minang: '+str(np.round(minang,0)), transform = ax.transAxes)
      ##ax.set_ylim(0,8)
      #ax.set_title('Triangle '+str(t))
      #ax.set_xlabel('ALS freeboard (m)')
      #ax.set_ylabel('Probability (%)')
      #pdfname = 'tri_'+str(t)
      #fig1.savefig(outpath_tri+pdfname)
      
      
      #select triangles where there was almost no area change
      #diff = (area2-area1)/1000000
      #if np.abs(diff) < .003:
      if np.abs(ddd*1e6) < .3:
	gtri1_nc.append(vert1)
	gtri2_nc.append(vert2)
	gtri1_nc_area.append(area1)
	gtri2_nc_area.append(area2)
	gtri1_nc_fb.extend(ttmp1.tolist())
	gtri2_nc_fb.extend(ttmp2.tolist())
	gtri_nc_div.append(ddd)
	
      else:
	gtri1_bc.append(vert1)
	gtri2_bc.append(vert2)
	gtri1_bc_area.append(area1)
	gtri2_bc_area.append(area2)
	#separate between increasing and decreasing area triangles
	if ddd>0:
	  gtri1_bc_pos.append(vert1)
	  gtri2_bc_pos.append(vert2)
	  gtri1_bc_pos_area.append(area1)
	  gtri2_bc_pos_area.append(area2)
	  mean1_bc_pos.append(np.mean(tmp1))
	  mean2_bc_pos.append(np.mean(tmp2))
	  gtri1_bc_pos_fb.extend(ttmp1.tolist())
	  gtri2_bc_pos_fb.extend(ttmp2.tolist())
	  gtri_bc_pos_div.append(ddd)
	  
	else:
	  gtri1_bc_neg.append(vert1)
	  gtri2_bc_neg.append(vert2)
	  gtri1_bc_neg_area.append(area1)
	  gtri2_bc_neg_area.append(area2)
	  mean1_bc_neg.append(np.mean(tmp1))
	  mean2_bc_neg.append(np.mean(tmp2))
	  gtri1_bc_neg_fb.extend(ttmp1.tolist())
	  gtri2_bc_neg_fb.extend(ttmp2.tolist())
	  gtri_bc_neg_div.append(ddd)

print 'done with the triangle business'

######################################################3
#Save stuff
np.save('gtri1',gtri1)
np.save('gtri2',gtri2)
np.save('gtri1_nc',gtri1_nc)
np.save('gtri2_nc',gtri2_nc)
np.save('gtri1_bc',gtri1_bc)
np.save('gtri2_bc',gtri2_bc)
np.save('gtri1_bc_pos',gtri1_bc_pos)
np.save('gtri2_bc_pos',gtri2_bc_pos)
np.save('gtri1_bc_neg',gtri1_bc_neg)
np.save('gtri2_bc_neg',gtri2_bc_neg)

np.save('mean1',mean1)
np.save('mean2',mean2)
np.save('mode1',mode1)
np.save('mode2',mode2)
np.save('bias',bias)

np.save('mean1_bc_pos',mean1_bc_pos)
np.save('mean2_bc_pos',mean2_bc_pos)
np.save('mean1_bc_neg',mean1_bc_neg)
np.save('mean2_bc_neg',mean2_bc_neg)

np.save('gtri1_area',gtri1_area)
np.save('gtri2_area',gtri2_area)
np.save('gtri1_nc_area',gtri1_nc_area)
np.save('gtri2_nc_area',gtri2_nc_area)
np.save('gtri1_bc_area',gtri1_bc_area)
np.save('gtri2_bc_area',gtri2_bc_area)
np.save('gtri1_bc_pos_area',gtri1_bc_pos_area)
np.save('gtri2_bc_pos_area',gtri2_bc_pos_area)
np.save('gtri1_bc_neg_area',gtri1_bc_neg_area)
np.save('gtri2_bc_neg_area',gtri2_bc_neg_area)

np.save('gtri1_nc_fb',gtri1_nc_fb)
np.save('gtri2_nc_fb',gtri2_nc_fb)
np.save('gtri1_bc_pos_fb',gtri1_bc_pos_fb)
np.save('gtri2_bc_pos_fb',gtri2_bc_pos_fb)
np.save('gtri1_bc_neg_fb',gtri1_bc_neg_fb)
np.save('gtri2_bc_neg_fb',gtri2_bc_neg_fb)

np.save('gtri_nc_div',gtri_nc_div)
np.save('gtri_bc_pos_div',gtri_bc_pos_div)
np.save('gtri_bc_neg_div',gtri_bc_neg_div)

dux = np.array(dux)
duy = np.array(duy)
dvx = np.array(dvy)
dvy = np.array(dvy)

np.save('dux',dux)
np.save('duy',duy)
np.save('dvx',dvx)
np.save('dvy',dvy)
np.save('gtri_minang',gtri_minang)
##exit()
#############################################################################3
#Load stuff
gtri1 = np.load('gtri1.npy')
gtri2 = np.load('gtri2.npy')
gtri1_area = np.load('gtri1_area.npy')
gtri2_area = np.load('gtri2_area.npy')

##to get just the not changing triangles
#gtri1 = np.load('gtri1_nc.npy')
#gtri2 = np.load('gtri2_nc.npy')
#outname1 = 'map_19_04_nc'
#outname2 = 'map_24_04_nc'
#outals1 = '../ALS_20150419_nc'
#outals2 = '../ALS_20150424_nc'

###to get just the changing triangles
##gtri1 = np.load('gtri1_bc.npy')
##gtri2 = np.load('gtri2_bc.npy')
##outname1 = 'map_19_04_bc'
##outname2 = 'map_24_04_bc'
##outals1 = '../ALS_20150419_bc'
##outals2 = '../ALS_20150424_bc'
##positive changes
#gtri1 = np.load('gtri1_bc_pos.npy')
#gtri2 = np.load('gtri2_bc_pos.npy')
#outname1 = 'map_19_04_bc_pos'
#outname2 = 'map_24_04_bc_pos'
#outals1 = '../ALS_20150419_bc_pos'
#outals2 = '../ALS_20150424_bc_pos'
##negative changes
#gtri1 = np.load('gtri1_bc_neg.npy')
#gtri2 = np.load('gtri2_bc_neg.npy')
#outname1 = 'map_19_04_bc_neg'
#outname2 = 'map_24_04_bc_neg'
#outals1 = '../ALS_20150419_bc_neg'
#outals2 = '../ALS_20150424_bc_neg'

dux = np.load('dux.npy')
duy = np.load('duy.npy')
dvx = np.load('dvx.npy')
dvy = np.load('dvy.npy')
mean1 = np.load('mean1.npy')
mean2 = np.load('mean2.npy')

#calculate deformation
div = dux + dvy
shr = .5*np.sqrt((dux-dvy)**2+(duy+dvx)**2)
dr = np.sqrt(div**2+shr**2)

np.save('div',div)
np.save('shr',shr)
np.save('dr',dr)

#Basemap map
fig1     = plt.figure(figsize=(10,10))
ax      = fig1.add_subplot(111)
ax.set_title(title1, fontsize=28)

m1 = pr.plot.area_def2basemap(area_def2)
m1.drawcoastlines()

#virtual buoys
ax.plot(x1,y1,'o',linewidth=2,color='purple')

mask = np.ones_like(als1, dtype=bool)
mask = False  
group = []
#separate triangles
for i in range(len(gtri1)):
    poly = shpol(gtri1[i])
    patch = PolygonPatch(poly, edgecolor='orchid', alpha=1, fill=False)
    ax.add_patch(patch)
    group.append(poly)
    
    #make mask
    px,py = poly.exterior.coords.xy
    pverts = np.vstack((px,py)).T
    path = Path(pverts)
    grid = path.contains_points(points1)
    grid = grid.reshape((xx1.shape[0],xx1.shape[1]))
    mask = np.logical_or(mask,grid)

###unified triangles   
#area1 = cascaded_union(group)
#patch = PolygonPatch(area1, edgecolor='k', alpha=1, fill=False)
#ax.add_patch(patch)


als_cut = np.ma.array(als1,mask=~mask)
#im1 = m1.pcolormesh(xx1, yy1, als_cut, cmap=plt.cm.jet, vmin=-0.2, vmax=0)
im1 = m1.pcolormesh(xx1, yy1, als_cut, cmap=plt.cm.jet, vmin=-0.2, vmax=1.2)
plt.colorbar(im1,orientation='horizontal')
#store the data
als_cut.dump(outals1)

# Then, I wanted to have a scale on the plot :
asb =  AnchoredSizeBar(ax.transData,
                          1000., # length of the bar in the data reference
                          "1 km", # label of the bar
                          loc=2, # location (lower right)
                          pad=0.1, borderpad=0.25, sep=5,
                          frameon=False)
ax.add_artist(asb) # and finaly, we add the scale to the inset axis


fig1.tight_layout()
fig1.savefig(outpath+outname1)
#exit()

fig2    = plt.figure(figsize=(10,10))
ax      = fig2.add_subplot(111)
ax.set_title(title2, fontsize=28)

m2 = pr.plot.area_def2basemap(area_def2)
m2.drawcoastlines()

#virtual buoys
ax.plot(x2,y2,'o',linewidth=2,color='purple')

mask = np.ones_like(als2, dtype=bool)
mask = False  
group = []
#separate triangles
for i in range(len(gtri2)):
    poly = shpol(gtri2[i])
    patch = PolygonPatch(poly, edgecolor='orchid', alpha=1, fill=False)
    ax.add_patch(patch)
    group.append(poly)
    
    #mask
    px,py = poly.exterior.coords.xy
    pverts = np.vstack((px,py)).T
    path = Path(pverts)
    grid = path.contains_points(points2)
    grid = grid.reshape((xx2.shape[0],xx2.shape[1]))
    mask = np.logical_or(mask,grid)
    
als_cut = np.ma.array(als2,mask=~mask)
#im2 = m2.pcolormesh(xx2, yy2, als_cut, cmap=plt.cm.jet, vmin=-0.2, vmax=0)
im2 = m2.pcolormesh(xx2, yy2, als_cut, cmap=plt.cm.jet, vmin=-0.2, vmax=1.2)
plt.colorbar(im2,orientation='horizontal')

###unified triangles   
#area2 = cascaded_union(group)
#patch = PolygonPatch(area2, edgecolor='k', alpha=1, fill=False)
#ax.add_patch(patch)

#store the data
als_cut.dump(outals2)

# Then, I wanted to have a scale on the plot :
asb =  AnchoredSizeBar(ax.transData,
                          1000., # length of the bar in the data reference
                          "1 km", # label of the bar
                          loc=2, # location (lower right)
                          pad=0.1, borderpad=0.25, sep=5,
                          frameon=False)
ax.add_artist(asb) # and finaly, we add the scale to the inset axis


fig2.tight_layout()
fig2.savefig(outpath+outname2)

##area change plot
#fig3    = plt.figure(figsize=(10,10))
#cx      = fig3.add_subplot(111)
#cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
#cx.text(.05, .95, 'c', ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

#m2 = pr.plot.area_def2basemap(area_def2)

##virtual buoys
#cx.plot(x2,y2,'o',linewidth=2,color='purple')

##triangles
#patches = []
#for i in range(len(gtri2)):
    #patch = Polygon(gtri2[i], edgecolor='orchid', alpha=1, fill=False)
    #patches.append(patch)
    
#diff = (np.array(gtri2_area) - np.array(gtri1_area))/1000000	#convert to km^2
#p = PatchCollection(patches, cmap=plt.cm.bwr, alpha=0.4)
#p.set_array(np.array(diff))
#cx.add_collection(p)
#p.set_clim([-.08, .08])
## create an axes on the right side of ax. The width of cax will be 5%
## of ax and the padding between cax and ax will be fixed at 0.05 inch.
#divider = make_axes_locatable(cx)
#cax = divider.append_axes("bottom", size="5%", pad=0.1)
#cbar = plt.colorbar(p, cax=cax,orientation='horizontal')
#cbar.set_label(r'Area change (km$^2$)',size=16)

#fig3.savefig(outpath+'map_area_change',bbox_inches='tight')

cmap=plt.cm.bwr
cdict = {'red':   ((0.0, 0.0, 0.0),
 		   (0.4, 1, 1),
		   (0.6, 1, 1),
		   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
		   (0.4, 1, 1),
		   (0.6, 1, 1),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 1.0),
		   (0.4, 1, 1),
		   (0.6, 1, 1),
                   (1.0, 0.0, 0.0))
        }

plt.register_cmap(name='MoreWhite', data=cdict)

#pick 1 out of 4!!!
deform = div*1e6
outname4 = 'map_div'
label = r'Divergence (10$^6$s$^{-1}$)'
interval = [-2.5, 2.5]
cmap=plt.cm.bwr
title = 'c'

#deform = shr*1e6
#outname4 = 'map_shr'
#label = r'Total shear (s$^{-1}$)'
#interval = [0, 3]
#cmap=plt.cm.Reds
#title = 'd'

#deform = dr*1e6
#outname4 = 'map_td'
#label = r'Deformation rate (s$^{-1}$)'
#interval = [0, 5]
#cmap=plt.cm.Reds
#title = 'e'

#rhow = 1025; rhoi = 917
#diff = np.array(mean2)-np.array(mean1)
#diff_t = rhow/(rhow-rhoi)*diff
#deform=diff
#outname4= 'map_diff_f'
#label = r'Freeboard difference (m)'
#interval = [-.2, .2]
#cmap = 'MoreWhite'
#cmap=plt.cm.bwr
#title = 'e'

#deform = bias
#outname4 = 'map_bias'
#label = r'Bias (m)'
#interval = [-.15, .15]
#cmap=plt.cm.bwr
#title = 'f'



#deformation plots
fig3    = plt.figure(figsize=(10,10))
cx      = fig3.add_subplot(111)
cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

m2 = pr.plot.area_def2basemap(area_def2)
m2.drawcoastlines()
#m2.drawparallels(np.arange(79.,90.,.01),labels=[1,0,0,0])
#m2.drawmeridians(np.arange(0.,360.,.1),latmax=90.,labels=[0,0,0,1,])

#virtual buoys
cx.plot(x2,y2,'o',linewidth=2,color='purple')

#triangles
print len(gtri2)
patches = []
for i in range(len(gtri2)):
    patch = Polygon(gtri2[i], edgecolor='orchid', alpha=1, fill=False)
    patches.append(patch)
    
p = PatchCollection(patches, cmap=cmap, alpha=0.4)
p.set_array(np.array(deform))
cx.add_collection(p)
p.set_clim(interval)
# create an axes on the right side of ax. The width of cax will be 5%
# of ax and the padding between cax and ax will be fixed at 0.05 inch.
divider = make_axes_locatable(cx)
cax = divider.append_axes("bottom", size="5%", pad=0.1)
cbar = plt.colorbar(p, cax=cax, orientation='horizontal')
cbar.set_label(label,size=16)

fig3.savefig(outpath+outname4,bbox_inches='tight')
##############################################################3
#outname5= 'map_diff_f'
#label = r'Freeboard difference (m)'
#interval = [-.3, .3]
#cmap=plt.cm.bwr

#cdict = {'red':   ((0.0, 0.0, 0.0),
 		   #(0.4, 1, 1),
		   #(0.6, 1, 1),
		   #(1.0, 1.0, 1.0)),

         #'green': ((0.0, 0.0, 0.0),
		   #(0.4, 1, 1),
		   #(0.6, 1, 1),
                   #(1.0, 0.0, 0.0)),

         #'blue':  ((0.0, 0.0, 1.0),
		   #(0.4, 1, 1),
		   #(0.6, 1, 1),
                   #(1.0, 0.0, 0.0))
        #}

#plt.register_cmap(name='MoreWhite', data=cdict)
#cmap = 'MoreWhite'

#rhow = 1025
#rhoi = 917
#diff = mean2-mean1#-0.05
#diff_t = rhow/(rhow-rhoi)*diff
##mean sea ice thickness change
#fig4    = plt.figure(figsize=(10,10))
#cx      = fig4.add_subplot(111)
#cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
#cx.text(.05, .95, title3, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

#m2 = pr.plot.area_def2basemap(area_def2)
#m2.drawcoastlines()
##m2.drawparallels(np.arange(79.,90.,.01),labels=[1,0,0,0])
##m2.drawmeridians(np.arange(0.,360.,.1),latmax=90.,labels=[0,0,0,1,])

##virtual buoys
#cx.plot(x2,y2,'o',linewidth=2,color='purple')

##triangles
#patches = []
#for i in range(len(gtri2)):
    #patch = Polygon(gtri2[i], edgecolor='orchid', alpha=1, fill=False)
    #patches.append(patch)
    
#p = PatchCollection(patches, cmap=cmap, alpha=.4)
#p.set_array(np.array(diff))
#cx.add_collection(p)
#p.set_clim(interval)
#cbar = plt.colorbar(p, orientation='horizontal')
#cbar.set_label(label,size=16)

#fig4.tight_layout()
#fig4.savefig(outpath+outname5,bbox_inches='tight')

exit()

###################################################################
#save shapefiles with both sets of buoys
#shapefiles with unified areas#####################################   
os.remove(outfile)
driver = ogr.GetDriverByName('Esri Shapefile')
ds = driver.CreateDataSource(outfile)
layer = ds.CreateLayer('', None, ogr.wkbPolygon)
# Add one attribute
layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
defn = layer.GetLayerDefn()

# Create a new feature (attribute and geometry)
feat = ogr.Feature(defn)
feat.SetField('id', 123)

# Make a geometry, from Shapely object
geom = ogr.CreateGeometryFromWkb(area.wkb)
feat.SetGeometry(geom)

layer.CreateFeature(feat)
feat = geom = None  # destroy these

#dalauney triangles shapefiles######################################
os.remove(outfile_tri)
driver = ogr.GetDriverByName('Esri Shapefile')
ds = driver.CreateDataSource(outfile_tri)
layer = ds.CreateLayer('', None, ogr.wkbPolygon)
# Add one attribute
layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
defn = layer.GetLayerDefn()

## If there are multiple geometries, put the "for" loop here
for i in group:
    # Create a new feature (attribute and geometry)
    feat = ogr.Feature(defn)
    feat.SetField('id', 123)

    # Make a geometry, from Shapely object
    geom = ogr.CreateGeometryFromWkb(i.wkb)
    feat.SetGeometry(geom)

    layer.CreateFeature(feat)
    feat = geom = None  # destroy these

# Save and close everything
ds = layer = feat = geom = None

#crop out the areas in commandline
#gdalwarp -cutline area_19_04.shp -crop_to_cutline -dstnodata "-9999.0" ../ALS_20150419_qgis.tiff ../ALS_20150419_cut.tiff
#gdalwarp -cutline area_24_04.shp -crop_to_cutline -dstnodata "-9999.0" ../ALS_20150424_qgis.tiff ../ALS_20150424_cut.tiff

