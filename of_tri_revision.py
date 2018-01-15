#! /usr/bin/python
# -*- coding: utf-8 -*-

#export virtual buoys vector layer created in QGIS to CSV files with geographical projection (WGS84, EPSG:4326)
#plot the ALS data and points

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
import osr, gdal, ogr
import pyresample as pr
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from scipy.interpolate import LinearNDInterpolator
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
outpath_tri = '../plots/pdf_tri_revision/'
#outpath_tri = '../plots/pdf_tri/'


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
vbuoys1 = '../vbuoys_1904_revision.csv'
vbuoys2 = '../vbuoys_2404_revision.csv'

vbuoys1 = '../vbuoys_1904_revision1.csv'
vbuoys2 = '../vbuoys_2404_revision1.csv'
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
#also mask the very high values (e.g. Lance, tents etc.)
als1 = np.ma.array(als1,mask=((als1==mv)))
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

als2 = np.ma.array(als2,mask=((als2==mv)))
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
trinum=[]
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

#less triangles (rev1)
offset[2]=.04
offset[4]=.08
offset[13]=.04
offset[12]=.04
offset[43]=.04
offset[96]=.02
offset[58]=-.1
offset[57]=-.08
offset[27]=-.06
offset[53]=-.05
offset[54]=-.07
offset[107]=-.04
offset[108]=-.02
offset[109]=-999
offset[60]=.0
offset[59]=-.04
offset[30]=-.06#visual interpolation between triangles 59 and 67 (along flight line) and conservation of mass
offset[67]=-.06
offset[66]=-.02
offset[68]=-.02
offset[15]=-.04#conservation of mass
offset[123]=-.04#conservation of mass
offset[69]=-.08
offset[65]=-.07#outlier
offset[31]=-.01#conservation of mass
offset[145]=-.03#conservation of mass
offset[144]=-.08#interpolation gives 0.03, but 'visual interpolation' along the flight line gives -.05, conservation of mass gives -.08
offset[140]=-.05
offset[55]=-.05
offset[122]=-.05
offset[131]=.0#check, conservation of mass
offset[133]=.01#conservation of mass
offset[132]=.02
offset[155]=-999#check, potential outlier!!!
offset[156]=.0
offset[153]=.0#conservation of mass
offset[89]=.02
offset[178]=.04
offset[179]=.02
offset[48]=.04
offset[49]=.04
offset[20]=.01
offset[21]=.01
offset[52]=.04
offset[22]=.0
offset[19]=.0
offset[26]=.02
offset[154]=.04
offset[175]=.04
offset[177]=.01
offset[181]=.03
offset[180]=.04
offset[182]=.09
offset[183]=.13
offset[185]=.09
offset[186]=.08
offset[211]=.1
offset[210]=.07
offset[207]=.08
offset[208]=.05
offset[209]=.07
offset[206]=.07
offset[202]=.05
offset[199]=.07
offset[203]=.02
offset[172]=.04
offset[119]=.01
offset[84]=.02
offset[85]=-.03
offset[40]=.02
offset[118]=.04
offset[200]=.05
offset[201]=.04
offset[151]=.02#check, edge value done by visual interpolation
offset[164]=.04
offset[163]=.02
offset[165]=.05
offset[168]=.02
offset[189]=.06
offset[187]=.05
offset[188]=.06
offset[190]=.04#conservation of mass
offset[191]=.01#conservation of mass - outlier, old value was 0.03
offset[192]=.02#conservation of mass
offset[171]=.01#conservation of mass
offset[161]=.01
offset[91]=.02
offset[82]=.02#interpolated value is high, beacuse its not along the flight line (close to 0.03)
offset[83]=.02#conservation of mass
offset[170]=.02 #iterp. along flight line
offset[169]=.03#conservation of mass - check! alternative: interpolate
offset[81]=.04
offset[146]=.07
offset[147]=.09#conservation of mass
offset[149]=.07#conservation of mass
offset[148]=.11#conservation of mass
offset[87]=.06
offset[88]=.07
offset[102]=.02#check, outlier, interpolate might be a better solution, but it will make it even bigger outlier#conservation of mass
offset[103]=.02#interpolated along flight line
offset[97]=.02#check, potential heavy outlier!!!
offset[98]=.03#check, interp. along flingt line
offset[99]=.06
offset[130]=.07
offset[129]=.1
offset[128]=.08
offset[167]=.07
offset[125]=.07
offset[127]=.11
offset[92]=.15
offset[72]=.1#conservation of mass
offset[73]=.12
offset[74]=.13
offset[94]=.08
offset[193]=.08
offset[194]=.09
offset[159]=.07
offset[160]=.05
offset[157]=.04#check
offset[158]=.05
offset[198]=.13
offset[197]=.11
offset[100]=.11
offset[101]=.12#check
offset[195]=.04
offset[93]=.02
offset[76]=.01
offset[120]=-.02#check
offset[79]=.01
offset[80]=.0
offset[113]=.01
offset[114]=.0
offset[115]=.01
offset[116]=.0
offset[117]=.01
offset[110]=.01
offset[136]=-999
offset[75]=.01#conservation of mass
offset[137]=.01#conservation of mass
offset[134]=.05#conservation of mass
offset[104]=.0#conservation of mass
offset[138]=.02#conservation of mass
offset[105]=.07
offset[106]=.08

##there is 0.5 more bias (regarding the negative freeboards applied) in the 1st overflights than in the 2nd
#offset = offset-.05

##new stuff
#offset[2]=.04
#offset[4]=.08
#offset[6]=-.05
#offset[12]=.04
#offset[13]=.04
#offset[14]=.04
#offset[20]=-.05
#offset[21]=-.1
#offset[27]=.0
#offset[28]=.01
#offset[29]=.01
#offset[30]=.0
#offset[34]=.02
#offset[35]=-.02
#offset[41]=.02
#offset[44]=.04
#offset[49]=-.02
#offset[50]=-.04
#offset[51]=-.08
#offset[53]=-.06
#offset[54]=-.07
#offset[55]=.04
#offset[56]=.04
#offset[59]=.04
#offset[60]=.02
#offset[62]=-.03
#offset[64]=-999
#offset[70]=.08
#offset[71]=.12
#offset[72]=.13
#offset[74]=.01
#offset[77]=.01
#offset[78]=.0
#offset[79]=.04
#offset[80]=.02#interpolated value is high, beacuse its not along the flight line (close to 0.03)
#offset[82]=.02
#offset[83]=-.03
#offset[85]=.06
#offset[86]=.07
#offset[87]=.02
#offset[88]=.15
#offset[89]=.02
#offset[90]=.08
#offset[91]=.02#check, potential heavy outlier!!!
#offset[92]=.03#check, potential heavy outlier!!!
#offset[93]=.06
#offset[94]=.11
#offset[95]=.12#check
#offset[96]=.03#check, outlier, interpolate might be a better solution, but it will make it even bigger outlier
#offset[97]=.03#interpolated value is high, beacuse its not along the flight line (close to 0.06)
#offset[99]=.07
#offset[100]=.08
#offset[103]=-.05
#offset[104]=-.05#interpolation gives 0.03, but 'visual interpolation' along the flight line gives -.05
#offset[105]=-999
#offset[106]=.0
#offset[107]=-.02
#offset[108]=-.13
#offset[110]=-.2
#offset[112]=-.05
#offset[114]=-.03
#offset[115]=.02
#offset[117]=.04
#offset[118]=.04
#offset[119]=.02
#offset[120]=.01
#offset[121]=.04
#offset[122]=-999#should be more negative, potential outlier!!!
#offset[123]=.0
#offset[124]=-.07
#offset[125]=-.02
#offset[126]=-.06
#offset[127]=-.02
#offset[128]=-.07
#offset[129]=.01
#offset[130]=-.1
#offset[131]=.01
#offset[132]=.01
#offset[133]=.0
#offset[134]=.01
#offset[135]=.0
#offset[136]=.01
#offset[137]=-.02
#offset[138]=.0
#offset[139]=.07
#offset[141]=.11
#offset[142]=.08
#offset[143]=.1
#offset[144]=.07
#offset[145]=-999
#offset[146]=.02
#offset[147]=-999
#offset[148]=.02
#offset[149]=.05
#offset[150]=.04
#offset[151]=.09
#offset[152]=.06
#offset[153]=.05
#offset[154]=.06
#offset[155]=.04
#offset[156]=.02
#offset[157]=.07
#offset[159]=.01
#offset[160]=.03
#offset[161]=.08
#offset[162]=.07
#offset[163]=-999
#offset[164]=.12
#offset[165]=.08
#offset[167]=.02#check, edge value done by visual interpolation
#offset[170]=.05
#offset[171]=.07
#offset[172]=.05
#offset[174]=.09
#offset[175]=.01
#offset[176]=.04
#offset[177]=.13
#offset[178]=.05
#offset[179]=.04
#offset[180]=.02
#offset[183]=.02
#offset[185]=.01
#offset[186]=-999#check
#offset[187]=-999#check
#offset[188]=-999#check
#offset[189]=-999#check
#offset[190]=-999
#offset[192]=.05
#offset[193]=.04
#offset[194]=.08
#offset[195]=.09
#offset[196]=.04
#offset[198]=.11
#offset[199]=.13
#offset[200]=.05
#offset[201]=.07
#offset[202]=.08
#offset[203]=.07
#offset[204]=.1
#offset[205]=.07
#offset[206]=.04
#offset[207]=.08
#offset[208]=.07
#offset[209]=.05
#offset[210]=.1
#offset[211]=.02
#offset[212]=.04
#offset[213]=.02




##################################################################################
#in this version we interpolate the bias values for the triangles that deformed

#first we have to read the deformation values form one of the old files (produced by the old version this routine of_tri.py towards the end)
#criteria has to be total deformation!!!
ddd = np.sqrt((np.load('dux.npy') + np.load('dvy.npy'))**2+(.5*np.sqrt((np.load('dux.npy')-np.load('dvx.npy'))**2+(np.load('duy.npy')+np.load('dvy.npy'))**2))**2)
mask = np.abs(ddd*1e6) > .3

#alternative option
#manually picked triangles that have nothing but thin ice and get deformed - hard to pin point the bias...
mask = np.load('bias.npy')==-999

#get the triangles and calculate centroids
centroids = np.mean(np.load('gtri2.npy'),axis=1)
#get the bias values
values = np.ma.compressed(np.ma.array(np.load('bias.npy'),mask=mask))

#mask deformed triangles
mask_tri = np.empty_like(centroids)
mask_triinv = np.empty_like(centroids)
for i in range(0,centroids.shape[1]):
    mask_tri[:,i]=mask
    mask_triinv[:,i]=~mask
centroids_nc = np.ma.compressed(np.ma.array(centroids,mask=mask_tri)).reshape(values.shape[0],2)

#interpolate just among the non-deformed triangles
g = LinearNDInterpolator(centroids_nc, values)
#get the values for the all the triangles
corbias = g(centroids)

#some values in the corners did not get interpolated (nan), use empiraical estimates for that
corbias = np.where(np.abs(corbias)>0,corbias,np.load('bias.npy'))

#alternative: interpolate hand picked triangles, where no thick ice exists

vt=-1
#print len(tripts1)
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
      trinum.append(t)
      
      #count valid traingles
      vt = vt+1
      
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
      #ttmp1 = ttmp1 + offset[t]
      ttmp1 = ttmp1 + corbias[vt]
      
      bias.append(offset[t])
      
      #get mean and mode
      me1 = np.mean(ttmp1)
      me2 = np.mean(ttmp2)
      mean1.append(me1)
      mean2.append(me2)
      
      #primary mode
      binval=200
      bins = np.linspace(-.2, 2,num=binval)
      count,lowerlimit,binsize,extra= stats.histogram(ttmp1,binval,defaultlimits=(-.1, 2))
      mindx1 =  np.argmax(count)
      count,lowerlimit,binsize,extra= stats.histogram(ttmp2,binval,defaultlimits=(-.1, 2))
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
      #ax.hist(ttmp1, bins, normed=True, histtype='stepfilled', alpha=0.4, lw=0, facecolor='red',  range=[-.1, 2])
      #ax.hist(ttmp2, bins, normed=True, histtype='stepfilled', alpha=0.4, lw=0, facecolor='blue',  range=[-.1, 2])
      #ax.text(.05,.95,'mean 1: '+str(me1), transform = ax.transAxes)
      #ax.text(.05,.9,'mean 2: '+str(me2), transform = ax.transAxes)
      #ax.text(.05,.85,'mean diff: '+str(me2-me1), transform = ax.transAxes)
      #ax.text(.05,.8,'div : '+str(np.round(ddd*1e6,2)), transform = ax.transAxes)
      #ax.text(.05,.75,'shr : '+str(np.round(sss*1e6,2)), transform = ax.transAxes)
      #ax.text(.05,.7,'offset: '+str(corbias[vt]), transform = ax.transAxes)
      #ax.text(.05,.65,'minang: '+str(np.round(minang,0)), transform = ax.transAxes)
      #me1_est = (me2*area2)/area1
      #ax.text(.05,.6,'mean 1 estimated: '+str(me1_est), transform = ax.transAxes)
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
np.save('trinum',trinum)
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
np.save('corbias',corbias)

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
vor = dvx-duy

np.save('div',div)
np.save('shr',shr)
np.save('dr',dr)

#Basemap map
fig1    = plt.figure(figsize=(10,10))
ax      = fig1.add_subplot(111)
ax.plot(.05, .95, 'w.', markersize=70, transform=ax.transAxes, markeredgecolor='k', markeredgewidth=2)
ax.text(.05, .95, title1, ha='center', va='center', transform=ax.transAxes, fontdict={'color':'k','size':30})

m1 = pr.plot.area_def2basemap(area_def1)
#virtual buoys
ax.plot(x1,y1,'o',linewidth=2,color='purple')
#Lance
lance_lon = 14.064979223402144;lance_lat = 83.122550903451426
xa, ya = m1(lance_lon, lance_lat)
ax.plot(xa,ya,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

mask = np.ones_like(als1, dtype=bool)
mask = False  
group = []
#separate triangles
for i in range(len(gtri1)):
    poly = shpol(gtri1[i])
    patch = PolygonPatch(poly, edgecolor='orchid', alpha=1, fill=False,linewidth=2)
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
m1.pcolormesh(xx1, yy1, als1, cmap=plt.cm.Greys, vmin=0, vmax=2)
im1 = m1.pcolormesh(xx1, yy1, als_cut, cmap=plt.cm.jet, vmin=0, vmax=1)
cbar = m1.colorbar(location='bottom')
cbar.set_label(r'Freeboard (m)',size=16)

#store the data
als_cut.dump(outals1)

#Scale
slon = 14.16;slat = 83.097
m1.drawmapscale(slon, slat, 10, 90, 1000, units='m', barstyle='fancy',fontsize=14)
#North arrow
xa, ya = m1(slon, slat)
ax.arrow(xa+800,ya,0,75,fc="k", ec="k", linewidth = 4, head_width=200, head_length=200,overhang=.5)
ax.text(xa+760,ya-150,'N')

fig1.savefig(outpath+outname1,bbox_inches='tight')

fig2    = plt.figure(figsize=(10,10))
ax      = fig2.add_subplot(111)
ax.plot(.05, .95, 'w.', markersize=70, transform=ax.transAxes, markeredgecolor='k', markeredgewidth=2)
ax.text(.05, .95, title2, ha='center', va='center', transform=ax.transAxes, fontdict={'color':'k','size':30})

m2 = pr.plot.area_def2basemap(area_def2)

#virtual buoys
ax.plot(x2,y2,'o',linewidth=2,color='purple')
#Lance
lance_lon = 15.331769876266121;lance_lat = 82.743808335535192
xa, ya = m2(lance_lon, lance_lat)
ax.plot(xa,ya,'*',markeredgewidth=2,color='hotpink',markersize=20,markeredgecolor='k')

mask = np.ones_like(als2, dtype=bool)
mask = False  
group = []
#separate triangles
for i in range(len(gtri2)):
    poly = shpol(gtri2[i])
    patch = PolygonPatch(poly, edgecolor='orchid', alpha=1, fill=False,linewidth=2)
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
m2.pcolormesh(xx2, yy2, als2, cmap=plt.cm.Greys, vmin=0, vmax=2)
im2 = m2.pcolormesh(xx2, yy2, als_cut, cmap=plt.cm.jet, vmin=0, vmax=1)
cbar = m2.colorbar(location='bottom')
cbar.set_label(r'Freeboard (m)',size=16)

###unified triangles   
#area2 = cascaded_union(group)
#patch = PolygonPatch(area2, edgecolor='k', alpha=1, fill=False)
#ax.add_patch(patch)

#store the data
als_cut.dump(outals2)

#Scale
slon = 15.43;slat = 82.718
m2.drawmapscale(slon, slat, 10, 90, 1000, units='m', barstyle='fancy',fontsize=14)
#North arrow
xa, ya = m2(slon, slat)
ax.arrow(xa+800,ya,0,75,fc="k", ec="k", linewidth = 4, head_width=200, head_length=200,overhang=.5)
ax.text(xa+760,ya-150,'N')

fig2.savefig(outpath+outname2,bbox_inches='tight')


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

deform = shr*1e6
outname4 = 'map_shr'
label = r'Total shear (s$^{-1}$)'
interval = [0, 3]
cmap=plt.cm.Reds
title = 'd'

#deform = vor*1e6
#outname4 = 'map_vor'
#label = r'Vorticity (s$^{-1}$)'
#interval = [-2.5, 2.5]
#cmap=plt.cm.bwr
#title = 'a'

#rhow = 1025; rhoi = 917
diff = np.array(mean2)-np.array(mean1)
#diff_t = rhow/(rhow-rhoi)*diff
deform=diff
outname4= 'map_diff_f'
label = r'Freeboard difference (m)'
interval = [-.2, .2]
cmap = 'MoreWhite'
cmap=plt.cm.bwr
title = 'e'

#deform = corbias
#outname4 = 'map_corbias'
#label = r'Freeboard correction (m)'
#interval = [-.15, .15]
#cmap=plt.cm.bwr
#title = 'f'



#deformation plots
fig3    = plt.figure(figsize=(10,10))
cx      = fig3.add_subplot(111)
if title != 'f':
  cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
  cx.text(.05, .95, title, ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

m2 = pr.plot.area_def2basemap(area_def2)
m2.drawcoastlines()

#virtual buoys
cx.plot(x2,y2,'o',linewidth=2,color='purple')

#triangles
print len(gtri2)
patches = []
for i in range(len(gtri2)):
    patch = Polygon(gtri2[i], edgecolor='orchid', alpha=1, fill=False)
    patches.append(patch)
    #label individual triangles
    #centroid = np.mean(gtri2[i],axis=0)
    #cx.text(centroid[0], centroid[1],trinum[i])
    
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

