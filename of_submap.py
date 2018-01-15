#! /usr/bin/python

import numpy as np
from mpl_toolkits.basemap import Basemap
import  matplotlib.pyplot as plt

outname='map_overview_sub'
path_out = '../plots/'

fig = plt.figure(figsize=(8,8))
ax1      = fig.add_subplot(111)

map1 = Basemap(projection='npstere',boundinglat=69,lon_0=0,resolution='l',round=True, ax=ax1)
map1.drawmapboundary(fill_color='#9999FF')
map1.fillcontinents(color='#ddaa66',lake_color='#9999FF')
map1.drawcoastlines()

#draw frame for the close up map
from matplotlib.patches import Polygon
lats = [ 81., 81., 83.5, 83.5 ]
lons = [ 5, 25, 25, 5 ]
x,y = map1(lons, lats)
xy = zip(x,y)
poly = Polygon( xy, edgecolor='purple', alpha=1, facecolor='red', linewidth=2)
ax1.add_patch(poly)

fig.tight_layout()
fig.savefig(path_out+outname,bbox_inches='tight')