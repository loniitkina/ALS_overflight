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

path = '../data/'
outpath = '../plots/'

#######################################33
#total area
totarea1 = np.load('gtri1_area.npy')
totarea2 = np.load('gtri2_area.npy')

atmp1 = np.load('gtri1_nc_fb.npy')
atmp2 = np.load('gtri2_nc_fb.npy')
btmp1 = np.load('gtri1_bc_pos_fb.npy')
btmp2 = np.load('gtri2_bc_pos_fb.npy')
ctmp1 = np.load('gtri1_bc_neg_fb.npy')
ctmp2 = np.load('gtri2_bc_neg_fb.npy')

tmp1 = np.append(atmp1,btmp1)
all_als1 = np.append(tmp1,ctmp1)
tmp1 = np.mean(all_als1)
print 'mean_freeboard als1: ', tmp1

tmp2 = np.append(atmp2,btmp2)
tmp2 = np.mean(np.append(tmp2,ctmp2))
print 'mean_freeboard als2: ', tmp2


#no divergence
als1 = np.load('gtri1_nc_fb.npy')
als2 = np.load('gtri2_nc_fb.npy')
outname = 'hist_nc_ddd'
area1 = np.load('gtri1_nc_area.npy')
area2 = np.load('gtri2_nc_area.npy')
div = np.mean(np.load('gtri_nc_div.npy'))*1e6
title = 'd'

##divergence
#als1 = np.load('gtri1_bc_pos_fb.npy')
#als2 = np.load('gtri2_bc_pos_fb.npy')
#outname = 'hist_bc_pos'
#area1 = np.load('gtri1_bc_pos_area.npy')
#area2 = np.load('gtri2_bc_pos_area.npy')
#div = np.mean(np.load('gtri_bc_pos_div.npy'))*1e6
#title = 'f'

#convergence
als1 = np.load('gtri1_bc_neg_fb.npy')
als2 = np.load('gtri2_bc_neg_fb.npy')
outname = 'hist_bc_neg'
area1 = np.load('gtri1_bc_neg_area.npy')
area2 = np.load('gtri2_bc_neg_area.npy')
div = np.mean(np.load('gtri_bc_neg_div.npy'))*1e6
title = 'e'

outpath = '../plots/'
########################################33
##read the ALS data
#als1 = np.concatenate(als1).ravel()
#als2 = np.concatenate(als2).ravel()

#print np.min(als1),np.max(als1)
#print np.min(als2),np.max(als2)

fb1 = np.mean(als1)
fb2 = np.mean(als2)
diff = fb2-fb1

#area
areasumtot = np.sum(totarea1)/1000000
print 'total area:',areasumtot
print 'area difference: ', np.sum(totarea2)/1000000 - np.sum(totarea1)/1000000

#subregion
areasum = np.sum(area1)/1000000
area_diff = (np.sum(area2)-np.sum(area1))/1000000
frac_diff = (area_diff/areasumtot)*100
frac_diff_loc = (area_diff/areasum)*100
print 'subregion area:',areasum
print 'area difference: ', area_diff

rhow = 1025
#rhoi = 895	#SYI at overflight site (gas chambers)
#rhoi = 917
rhoi = 906	#130 cm, partially de-brined
it_diff = (rhow/(rhow-rhoi))*(fb2-fb1)
it_diff_all = (rhow/(rhow-rhoi))*(tmp2-tmp1)
print 'ice thickness change (m): ', it_diff
print 'ice thickness change (m) - whole area: ', it_diff_all
#exit()

#get the mode
from scipy.stats import mstats
mode_als1 = mstats.mode(als1, axis=None)[0][0]
mode_als2 = mstats.mode(als2, axis=None)[0][0]
print 'mode ALS1:', mode_als1
print 'mode ALS2:', mode_als2

#get thin ice mode
als1_thin = np.ma.array(als1,mask=als1>.18)
als2_thin = np.ma.array(als2,mask=als2>.18)
mode_als1_thin = mstats.mode(als1_thin, axis=None)[0][0]
mode_als2_thin = mstats.mode(als2_thin, axis=None)[0][0]
print 'mode ALS1 thin ice:', mode_als1_thin
print 'mode ALS2 thin ice:', mode_als2_thin

fig1    = plt.figure(figsize=(8,8))
ax      = fig1.add_subplot(111)
ax.plot(.06, .94, 'w.', markersize=80, transform=ax.transAxes, markeredgecolor='k', markeredgewidth=2)
ax.text(.06, .94, title, ha='center', va='center', transform=ax.transAxes, fontdict={'color':'k','size':35})

binno = 171
x = als1#.compressed().flatten()
#print help(ax.hist)
n, bins, patches = ax.hist(x, binno, normed=True, histtype='step', color='purple', alpha=.8, range=[-.2, 1.5], label='19/04/2015', lw = 3)

x = als2#.compressed().flatten()
n, bins, patches = ax.hist(x, binno, normed=True, histtype='step', color='royalblue', alpha=.8, range=[-.2, 1.5], label='24/04/2015', lw = 3)

text = 'area difference: '+ str(round(area_diff,2)) + r' km$^2$'
ax.text(.43, .8, text, ha='left', va='center', transform=ax.transAxes, fontdict={'color':'k','size':16})
text = 'area fraction: ' + str(round(frac_diff,0)) + ' %' + ' (' + str(round(frac_diff_loc,0)) + '%)'
ax.text(.43, .75, text, ha='left', va='center', transform=ax.transAxes, fontdict={'color':'k','size':16})
text = r'$\bar{h}_f $ on 19th: '+ str(round(fb1,2)) + r' m'
ax.text(.43, .7, text, ha='left', va='center', transform=ax.transAxes, fontdict={'color':'k','size':16})
text = r'$\bar{h}_f $ on 24th: '+ str(round(fb2,2)) + r' m'
ax.text(.43, .65, text, ha='left', va='center', transform=ax.transAxes, fontdict={'color':'k','size':16})
text = r'$\Delta \bar{h} $: '+ str(round(it_diff,2)) + r' m'
ax.text(.43, .6, text, ha='left', va='center', transform=ax.transAxes, fontdict={'color':'k','size':16})
text = r'$\bar{\epsilon}_{div}$: '+ str(round(div,2)) + r' 10$^{-6}$s$^{-1}$'
ax.text(.43, .55, text, ha='left', va='center', transform=ax.transAxes, fontdict={'color':'k','size':16})

ax.set_xlim(-.2, 1.5)
ax.set_ylim(0, 3)
ax.set_xlabel('Freeboard (m)',size=18)
ax.set_ylabel('Probability (%)',size=18)
plt.xticks(size=16)
plt.yticks(size=16)
ax.legend(loc='upper right',prop={'size':16}, fancybox=True, framealpha=.5, numpoints=1)


#estimate deformed ice volume
if title == 'e':
  fig2    = plt.figure(figsize=(8,8))
  bx      = fig2.add_subplot(111)
  
  rhow = 1025
  #rhoi = 895	#SYI at overflight site (gas chambers)
  #rhoi = 917
  rhoi = 906	#130 cm, partially de-brined
  it1 = (rhow*als1)/(rhow-rhoi)
  it2 = (rhow*(als2+.0))/(rhow-rhoi)

  vol1 = it1*25	#5x5m pixel
  vol2 = it2*25

  #n1, bins, patches = plt.hist(vol1, binno, normed=False, label='19/04/2015', histtype='step', range=[20, 250])
  #n2, bins, patches = plt.hist(vol2, binno, normed=False, label='24/04/2015', histtype='step', range=[20, 250])
  
  n1, bins, patches = plt.hist(als1, binno, normed=False, label='19/04/2015', histtype='step', range=[.1, 1.5])
  n2, bins, patches = plt.hist(als2, binno, normed=False, label='24/04/2015', histtype='step', range=[.1, 1.5])
  n3, bins, patches = plt.hist(all_als1, binno, normed=False, label='19/04/2015', histtype='step', range=[.1, 1.5])
  
  bx.legend(loc='upper right',prop={'size':16})
  #bx.set_ylim(-2,2)
  x = (bins[1:] + bins[:-1])/2
  print n1
  print n2
  bx.plot(x,n2-n1)
  fig2.savefig(outpath+'hist_deform',bbox_inches='tight')
  
  voldiff = np.sum(x*n2) - np.sum(x*n1)
  print 'freeboard difference in deformed ice (m): '+str(voldiff)
  print 'freeboard/thickness/volume increase of deformed ice (%): '+str(voldiff/np.sum(x*n1)*100)
  print 'freeboard/thickness/volume increase of deformed ice in total volume of thick ice (%): '+str(voldiff/np.sum(x*n3)*100)
  #winter 4 storms/month for 9 months (Oct-Jun) = 36 storms = 11.7% as a lower bound (no major storm), FS=30%
  text = 'def. ice volume change: '+ str(round(voldiff/np.sum(x*n3)*100,2)) + ' %'
  ax.text(.43, .5, text, ha='left', va='center', transform=ax.transAxes, fontdict={'color':'k','size':16})
fig1.savefig(outpath+outname,bbox_inches='tight')