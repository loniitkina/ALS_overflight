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

title = 'ICE-ARC overflight 19-24/04/2015'

div = np.load('div.npy')*1e6#*3600*24*30.5	#convert from s^-1 to mon^-1
shr = np.load('shr.npy')*1e6#*3600*24*30.5	#convert from s^-1 to mon^-1
dr = np.load('dr.npy')*1e6#*3600*24*30.5	#convert from s^-1 to mon^-1
mean1 = np.load('mean1.npy')
mean2 = np.load('mean2.npy')
#mode1 = np.load('mode1.npy')
#mode2 = np.load('mode2.npy')
#mode1b = np.load('mode1b.npy')
#mode2b = np.load('mode2b.npy')
area1 = np.load('gtri1_area.npy')
area2 = np.load('gtri2_area.npy')
minang = np.load('gtri_minang.npy')

print 'total area1: ', np.sum(area1)/1000**2
print 'total area2: ', np.sum(area2)/1000**2
print 'triangle sizes: min, avg, max: ', np.min(area1),np.mean(area1), np.max(area1)
#volume difference
#vol1 = area1*np.ma.array(mean1,mask=minang<25)
#vol1 = area1*np.ma.array(mean1,mask=area2<50000)

rhow = 1025
#rhoi = 895	#SYI at overflight site (gas chambers)
#rhoi = 917
rhoi = 906	#130 cm, partially de-brined
diff_all = (mean2-mean1)#/mean1*100
diff_a = (area2-area1)
it1 = (rhow/(rhow-rhoi))*(mean1)
it2 = (rhow/(rhow-rhoi))*(mean2)
it_diff = (rhow/(rhow-rhoi))*(mean2-mean1)

diff_vol = (it2*area2)-(it1*area1)
vol1 = area1*it1
vol2 = area2*it2
#print 'volume difference: ', np.sum(vol2-vol1), 'm^3'
print 'volume difference: ', np.sum(diff_vol), 'm^3'
#print 'relative increase: ', np.round(((vol2-vol1)/vol1)*100,0)
#print 'relative increase, area: ', np.round(((area2-area1)/area1)*100,0)
print 'relative increase, volume: ', np.sum(vol2-vol1)/np.sum(vol1), '%'
print 'relative increase, area: ', np.sum(area2-area1)/np.sum(area1), '%'
print 'increase, fb: ', np.sum(mean2-mean1), 'm'
print 'increase, thickness: ', np.sum(it_diff), 'm'
print 'increase, ice volume: ', (np.sum(it_diff)*np.sum(area1))/1e9, 'km^3'
print 'realtive increase, ice volume (mean it=2m): ', (np.sum(it_diff)*np.sum(area1))/(np.sum(area1)*2)

#plt.bar(range(0,len(vol1)),np.round(((vol2-vol1)/vol1)*100,0),label='volume')
#plt.plot(range(0,len(vol1)),np.round(((area2-area1)/area1)*100,0),'o',label='area',color='g')
#plt.plot(range(0,len(vol1)),np.round(((mean2-mean1)/mean1)*100,0),'o',label='freeboard',color='r')
#plt.ylim(-50,100)
#plt.legend()
#plt.show()

#plt.plot(range(0,len(vol1)),area1,label='area1')
#plt.plot(range(0,len(vol1)),area2,label='area2')
#plt.legend()
#plt.show()


#plt.plot(range(0,len(vol1)),diff_vol,label='volume difference')
##plt.plot(range(0,len(vol1)),mean2-mean1,label='freeboard difference')
#plt.legend()
#plt.show()







#select just those FB where the triangles are large (after deformation)
mask = (area2<5000)
#get just those that satifsy the minang criteria
mask = (minang<1)
#get just those where freeboard change is connected to the deformation

diff = np.ma.array(diff_all,mask=mask)
shr = np.ma.array(shr,mask=mask)
div = np.ma.array(div,mask=mask)
diff_it = np.ma.array(it_diff,mask=mask)
mean1 = np.ma.array(mean1,mask=mask)
diff = np.ma.compressed(diff)
shr = np.ma.compressed(shr)
div = np.ma.compressed(div)
diff_it = np.ma.compressed(diff_it)
mean1 = np.ma.compressed(mean1)

#divergence
fig1    = plt.figure(figsize=(8,11))
ax      = fig1.add_subplot(111)
ax.plot(.05, .95, 'w.', markersize=70, transform=ax.transAxes, markeredgecolor='k', markeredgewidth=2)
ax.text(.05, .95, 'a', ha='center', va='center', transform=ax.transAxes, fontdict={'color':'k','size':30})
ax.set_xlabel(r'$\epsilon_{div} (10^6s^{-1}$)', fontsize=20)
ax.set_ylabel(r'$\Delta h_f (m)$', fontsize=20)
#ax.scatter(div,diff_all,c='.8')
#ax.plot(div,diff, 'o',markerfacecolor='skyblue',ms=8,mew=1)



ax.set_ylim(-.3,.3)
ax.set_xlim(-2.5,2.5)
sc = ax.scatter(div,diff,c=mean1,marker='o',s=80,alpha=0.8, cmap=plt.cm.bwr)
cbar = plt.colorbar(sc, orientation='horizontal')
cbar.ax.set_xlabel(r'$h_f (m)$',size=20)

##only non-shearing cases
#div_ns = np.ma.array(div,mask=shr>.8)
#ax.scatter(div_ns,diff,marker='x',s=80, c='k')



#divider = make_axes_locatable(ax)
#cax = divider.append_axes("b", size="5%", pad=0.1)
#cbar = plt.colorbar(p, cax=cax, orientation='horizontal')
#cbar.set_label(label,size=16)



#divergence and sea ice thickness
#just plot axis, adjust to the min,max and dont plot any data (avoid overlaying)
aax = ax.twinx()
aax.spines['top'].set_visible(False)
aax.spines['bottom'].set_visible(False)
aax.set_ylabel(r'$\Delta h (m)$', fontsize=20)
amin,amax=ax.get_ylim()
aamin = (rhow/(rhow-rhoi))*amin
aamax = (rhow/(rhow-rhoi))*amax
aax.set_ylim(aamin,aamax)
aminx,amaxx=ax.get_xlim()
aax.set_xlim(aminx,amaxx)
#aax.plot(div,diff_it, 'o',markerfacecolor='skyblue',ms=8,mew=1)




#bin the data
from scipy import stats
binval = np.arange(-2,2.5,.5)
bin_means, bin_edges, binnumber = stats.binned_statistic(div, diff, statistic='mean', bins=binval)
x = (bin_edges[1:] + bin_edges[:-1])/2

fit = np.polyfit(x,bin_means,1)
fit = np.polyfit(div,diff,1)
fit_fn = np.poly1d(fit) 

#print fit
#y = fit[0]*x+fit[1] # regression line
#print x
#print bin_means
ax.plot(x, fit_fn(x), '--k')#,label='best fit')
#ax.plot(x,bin_means,'*',markerfacecolor='r',ms=20,mew=1,label='means of bins')


##Pearson R
#print np.corrcoef(div,diff)[1,0]

from scipy.stats import pearsonr
[pr,p] = pearsonr(div,diff)
r2 = r'$R= %1.2f$' %(pr**2)
linreg = r'$\Delta h_f = %1.3f*\epsilon_{div} %1.3f$' %(fit[0],fit[1])
ax.text(.05, .1, r2, ha='left', va='center', transform=ax.transAxes, fontdict={'color':'k','size':18})
ax.text(.05, .05, linreg, ha='left', va='center', transform=ax.transAxes, fontdict={'color':'k','size':18})



#xi = np.array([ x, np.ones(len(x))])
#w = np.linalg.lstsq(xi.T,bin_means)[0] # obtaining the parameters
#print w
## plotting the line
#line = w[0]*xi+w[1] # regression line
#ax.plot(xi,line,'r-')

#exit()


ax.legend(loc='upper right',prop={'size':18}, fancybox=False, framealpha=0.5,numpoints=1)
fig1.savefig('../plots/corr_div',bbox_inches='tight')

#shear
fig2    = plt.figure(figsize=(8,8))
bx      = fig2.add_subplot(111)
bx.plot(.05, .95, 'w.', markersize=70, transform=bx.transAxes, markeredgecolor='k', markeredgewidth=2)
bx.text(.05, .95, 'b', ha='center', va='center', transform=bx.transAxes, fontdict={'color':'k','size':30})
bx.set_xlabel(r'Shear (10$^6$s$^{-1}$)', fontsize=24)
bx.set_ylabel('Mean freeboard change (m)', fontsize=24)
#bx.scatter(shr,diff_all,c='.8')
#bx.plot(shr,diff, 'bo', shr, fit_fn(shr), '--k',markerfacecolor='skyblue')
#bx.plot(shr,diff, 'o',markerfacecolor='skyblue',ms=8,mew=1)


sc = bx.scatter(shr,diff,c=mean1,marker='o',s=80,alpha=0.8, cmap=plt.cm.bwr)
cbar = plt.colorbar(sc)
cbar.ax.set_ylabel(r'$h_f (m)$',size=24)



#bin the data
from scipy import stats
binval = np.arange(-.25,2.5,.5)
bin_means, bin_edges, binnumber = stats.binned_statistic(shr, diff, statistic='mean', bins=binval)
x = (bin_edges[1:] + bin_edges[:-1])/2

fit = np.polyfit(x,bin_means,1)
fit_fn = np.poly1d(fit) 

print x
print bin_means
bx.plot(x,bin_means,'*',markerfacecolor='red',ms=20,mew=1,label='means of bins')
bx.plot(x, fit_fn(x), '--k',label='best fit')


[pr,p] = pearsonr(shr,diff)
r2 = r'$R= %1.2f$' %(pr**2)
linreg = r'$\Delta h_f = %1.3f*\epsilon_{shr} %1.3f$' %(fit[0],fit[1])
bx.text(.05, .1, r2, ha='left', va='center', transform=bx.transAxes, fontdict={'color':'k','size':20})
bx.text(.05, .05, linreg, ha='left', va='center', transform=bx.transAxes, fontdict={'color':'k','size':20})




fig2.savefig('../plots/corr_shr',bbox_inches='tight')

#multiple lienar regression
import statsmodels.api as sm
def reg_m(y, x):
    ones = np.ones(len(x[0]))
    X = sm.add_constant(np.column_stack((x[0], ones)))
    for ele in x[1:]:
        X = sm.add_constant(np.column_stack((ele, X)))
    results = sm.OLS(y, X).fit()
    return results

y = diff
x = [div,shr]
print reg_m(y, x).summary()


#just shear at convergence
shr = np.load('shr.npy')*1e6
div = np.load('div.npy')*1e6
mean1 = np.load('mean1.npy')
mask = (div>0) | (minang<1) #| (area2<40000)
diff  = np.ma.array(diff_all,mask=mask)
shr = np.ma.array(shr,mask=mask)
div = np.ma.array(div,mask=mask)
diff_it  = np.ma.array(it_diff,mask=mask)
mean1  = np.ma.array(mean1,mask=mask)
diff = np.ma.compressed(diff)
shr = np.ma.compressed(shr)
div = np.ma.compressed(div)
diff_it = np.ma.compressed(diff_it)
mean1 = np.ma.compressed(mean1)

fig3    = plt.figure(figsize=(8,11))
cx      = fig3.add_subplot(111)
cx.plot(.05, .95, 'w.', markersize=70, transform=cx.transAxes, markeredgecolor='k', markeredgewidth=2)
cx.text(.05, .95, 'b', ha='center', va='center', transform=cx.transAxes, fontdict={'color':'k','size':30})

cx.set_xlabel(r'$\epsilon_{shr} (10^6s^{-1}$)', fontsize=20)
cx.set_ylabel(r'$\Delta h_f (m)$', fontsize=20)
#bx.scatter(dr,diff_all,c='.8')
#bx.plot(dr,diff, 'bo', dr, fit_fn(dr), '--k',markerfacecolor='skyblue')
#cx.plot(shr,diff, 'o',markerfacecolor='skyblue',ms=8,mew=1)

cx.set_ylim(-.02,.3)
cx.set_xlim(-.05,3.05)
sc = cx.scatter(shr,diff,c=mean1,marker='o',s=80,alpha=0.8, cmap=plt.cm.bwr)
cbar = plt.colorbar(sc, orientation='horizontal')
cbar.ax.set_xlabel(r'$h_f (m)$',size=20)





#shear and sea ice thickness
#just plot axis, adjust to the min,max and dont plot any data (avoid overlaying)
ccx = cx.twinx()
ccx.spines['top'].set_visible(False)
ccx.spines['bottom'].set_visible(False)
ccx.set_ylabel(r'$\Delta h (m)$', fontsize=20)
cmin,cmax=cx.get_ylim()
ccmin = (rhow/(rhow-rhoi))*cmin
ccmax = (rhow/(rhow-rhoi))*cmax
ccx.set_ylim(ccmin,ccmax)
cminx,cmaxx=cx.get_xlim()
ccx.set_xlim(cminx,cmaxx)
#ccx.plot(shr,diff_it, 'o',markerfacecolor='skyblue',ms=8,mew=1)




#bin the data
from scipy import stats
binval = np.arange(-.25,3,.5)
bin_means, bin_edges, binnumber = stats.binned_statistic(shr, diff, statistic='mean', bins=binval)
x = (bin_edges[1:] + bin_edges[:-1])/2

fit = np.polyfit(shr,diff,1)
fit_fn = np.poly1d(fit) 


[pr,p] = pearsonr(shr,diff)
#print p
#exit()
r2 = r'$R= %1.2f$' %(pr**2)
linreg = r'$\Delta h_f = %1.3f*\epsilon_{shr}+%1.3f$' %(fit[0],fit[1])
cx.text(.5, .1, r2, ha='left', va='center', transform=cx.transAxes, fontdict={'color':'k','size':16})
cx.text(.5, .05, linreg, ha='left', va='center', transform=cx.transAxes, fontdict={'color':'k','size':16})


print x
print bin_means
cx.plot(x, fit_fn(x), '--k')#,label='best fit')
#cx.plot(x,bin_means,'*',markerfacecolor='red',ms=20,mew=1,label='means of bins')
fig3.savefig('../plots/corr_shr_conv',bbox_inches='tight')


#multiple linear regression
y = diff
x = [div,shr]
print reg_m(y, x).summary()