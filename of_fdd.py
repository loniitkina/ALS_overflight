#! /usr/bin/python

import numpy as np
from datetime import datetime
from of_func import *

#ALS overflights
start = datetime(2015, 4, 19, 0, 0)
end = datetime(2015, 4, 24, 0, 0)

#winter
start = datetime(2015, 1, 15, 0, 0)
end = datetime(2015, 4, 24, 0, 0)



#meteorological data (from Lance's met system)
metfile = '../data/10minute_nounits.csv'
mettime = getColumn(metfile,0)
metdates = [ datetime.strptime(mettime[x], "%Y-%m-%d %H:%M:%S") for x in range(len(mettime)) ]
mettemp = np.asarray(getColumn(metfile,6),dtype=float)
mettemp = np.ma.masked_invalid(mettemp)

#average the wind speed for every day
n = 6*24
col = len(mettemp)/n
tmp = mettemp.reshape(col,n)
atmp = np.mean(tmp,axis=1)

#get the hourly dates
dtmp = metdates[::n]
wsi = np.argmin(abs(np.asarray(dtmp)-start))
wei = np.argmin(abs(np.asarray(dtmp)-end))

wtemp = atmp[wsi:wei]
print 'mean winter temperature:'
print np.mean(wtemp)

##constant temp
#wtemp=np.ones_like(wtemp)*-20
##random 2 weeks
#wtemp=np.ones(7)*-20
#print wtemp

t0= -1.9

fdd = np.sum(t0-wtemp)

#FDD model (results in cm), (Lebedev, 1938)
#no initial thickness
ic = 1.33 * fdd ** 0.58
print ic

#initial thickness
it=0
ic = np.sqrt(it**2+(fdd/804.2082)) 
print ic
#it=5
#ic = np.sqrt(it**2+(fdd/804.2082)) 
#print ic-it
#it=20
#ic = np.sqrt(it**2+(fdd/804.2082)) 
#print ic-it
#it=50
#ic = np.sqrt(it**2+(fdd/804.2082)) 
#print ic-it
#it=100
#ic = np.sqrt(it**2+(fdd/804.2082)) 
#print ic-it