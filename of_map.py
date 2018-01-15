#! /usr/bin/python

#map the full lance track (or just one lance track)

import numpy as np
from mpl_toolkits.basemap import Basemap
import  matplotlib.pyplot as plt
from datetime import datetime, timedelta
from glob import glob


dep='2'
letter = 'a) '

lance = '../data/lance_f3_of.csv'
als_track = '../flight_tracks'

b1 = '../buoys_1904.csv'
b2 = '../buoys_2404.csv'
outname='map_overview'
path_out = '../plots/'


fig = plt.figure(figsize=(12,12))
ax      = fig.add_subplot(111)
#ax.plot(.03, .95, 'w.', markersize=90, transform=ax.transAxes, markeredgecolor='k', markeredgewidth=2)
#ax.text(.03, .95, 'a', ha='center', va='center', transform=ax.transAxes, fontdict={'color':'k','size':40})

#ax1 = plt.subplot2grid((2,3), (1,0))
#plt = plt.subplot2grid((2,3), (0,1), colspan=2, rowspan=2)



#map1 = Basemap(projection='npstere',boundinglat=69,lon_0=0,resolution='l',round=True, ax=ax1)
#map1.drawmapboundary(fill_color='#9999FF')
#map1.fillcontinents(color='#ddaa66',lake_color='#9999FF')
#map1.drawcoastlines()

##draw frame for the close up map
#from matplotlib.patches import Polygon
#lats = [ 81., 81., 83.5, 83.5 ]
#lons = [ 5, 25, 25, 5 ]
#x,y = map1(lons, lats)
#xy = zip(x,y)
#poly = Polygon( xy, edgecolor='purple', alpha=1, facecolor='red', linewidth=2)
#ax1.add_patch(poly)


m = Basemap(resolution='h',
	    projection='stere', lat_ts=82, lat_0=90., lon_0=0.,
	    llcrnrlon=9, llcrnrlat=82.6,
	    urcrnrlon=20, urcrnrlat=83.05)

m.drawmapboundary(fill_color='#9999FF')
##Draw parallels and meridians
m.drawparallels(np.arange(80.,90.,.25),labels=[1,0,0,0], fontsize=16)
m.drawmeridians(np.arange(0.,30.,2.),latmax=85.,labels=[0,0,0,1], fontsize=16)
m.drawmapscale(17.2, 82.4, 15, 83, 20, barstyle='fancy', units='km', fontsize=16, yoffset=None, labelstyle='simple', fontcolor='k', fillcolor1='w', fillcolor2='k', ax=None, format='%d', zorder=None)

##plot ALS flight lines
als = m.readshapefile(als_track, 'als')
cnt = 0
for info, xy in zip(m.als_info, m.als):
    if info['flight'] == 'LYR-Lance-LYR':
      color = 'darkblue'
      label='ALS overflight 19/04/2015'
    else:
      color = 'k'
      label='ALS overflight 24/04/2015'
      cnt = cnt+1
    if cnt == 0:
      m.plot(xy[0], xy[1], lw=3, marker='o', label=label, markeredgecolor=color, color=color); cnt=cnt+1
    elif cnt == 3:
      m.plot(xy[0], xy[1], lw=3, marker='o', label=label, markeredgecolor=color, color=color)
    else:
      m.plot(xy[0], xy[1], lw=3, marker='o', markeredgecolor=color, color=color)

#plot Lance track
data = np.loadtxt(lance, dtype=np.str,delimiter=',', skiprows=1)
latitude = np.ma.array(data[:,1],mask=data[:,1]=='0')
longitude = np.ma.array(data[:,2],mask=data[:,2]=='0')
date = data[:,0]
x, y = m(longitude, latitude)
ax.plot(x,y,'-',linewidth=4, label='R/V Lance drift track',color='purple')

#Filter out the milestone dates and plot as additional points    
for i in range(0,len(x)):
  d = datetime.strptime(date[i], "%Y-%m-%d %H:%M:%S")
  event = datetime(2015,4,19,12,0)
  if event - timedelta(seconds=30*60) < d < event + timedelta(seconds=30*60):
    ax.plot(x[i],y[i],'o',color='goldenrod',markersize=18,markeredgewidth=1,label='R/V Lance 19/04/2015')
  event = datetime(2015,4,24,12,0)
  if event - timedelta(seconds=30*60) < d < event + timedelta(seconds=30*60):
    ax.plot(x[i],y[i],'o',color='r',markersize=18,markeredgewidth=1,label='R/V Lance 24/04/2015')

#Plot buoy positons
#data = np.loadtxt(b1, dtype=np.str,delimiter=',', skiprows=1)
#latitude = np.ma.array(data[:,2])
#longitude = np.ma.array(data[:,3])
#x, y = m(longitude, latitude)
#ax.plot(x,y,'o',markersize=10,markeredgewidth=1, label='GPS drifters',color='orchid')

cnt = 0
data = np.loadtxt(b2, dtype=np.str,delimiter=',', skiprows=1)
latitude = np.ma.array(data[:,2])
longitude = np.ma.array(data[:,3])
x, y = m(longitude, latitude)
if cnt==0:
  ax.plot(x,y,'o',markersize=10,markeredgewidth=1,color='hotpink',label='GPS drifters 24/04/2015'); cnt=cnt+1
else:
  ax.plot(x,y,'o',markersize=10,markeredgewidth=1,color='hotpink')


lg = ax.legend(loc='upper right',prop={'size':16}, fancybox=True, framealpha=.5, numpoints=1)
lg.get_frame().set_facecolor('white')

fig.tight_layout()
fig.savefig(path_out+outname,bbox_inches='tight')