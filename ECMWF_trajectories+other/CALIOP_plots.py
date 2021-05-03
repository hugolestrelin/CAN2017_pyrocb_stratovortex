#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Hugo Lestrelin

"""
from __future__ import print_function
import sys,os  
import socket      
if 'ciclad' in socket.gethostname():
    import sys,os
    sys.path.append('/home/hlestrel/pylib')                
from datetime import datetime,timedelta
from ECMWF_N import ECMWF 
import numpy as np 
import matplotlib.pyplot as plt 
import pickle,gzip
import matplotlib.colors as colors
import glob
from astropy.time import Time, TimeDelta
import matplotlib.colors as colors
import geosat
from cartopy import feature
import cartopy.crs as ccrs
import argparse
import scipy.signal as ss

## fonction transform proj -> legend
def proj2legend(x,y,x0,x1,y0,y1):
    xproj=x1-x0
    yproj=y1-y0
    xlegend=x/xproj
    ylegend=y/yproj
    return(xlegend,ylegend)

initials='CAN2017'
savfile='/data/hlestrel/fig_bulle_'+initials+'/track_CALIOP/'
Qs = 5.167*1.e-31
kbw = 1.0313
medfilter=True
widthfilter = 21
vmax = 15
vmin=1.5
lvlcon=np.arange(-20,-4.5,15)
cmap = 'gist_ncar'
# load dictionary of traces
with gzip.open(savfile+'selCaliop_'+initials+'.pkl','rb') as f:
    box,idate0,nbday,Cald = pickle.load(f)
with gzip.open(savfile+'export','rb') as f:
    export=pickle.load(f)
with gzip.open(savfile+'Ntrack-L1-'+initials,'rb') as f:
    CALIOPtrack=pickle.load(f)
Ctracklon=[]
Ctracklat=[]
Ctrackalt=[]
for i in range(len(CALIOPtrack)):
    Ctrackalt.append(CALIOPtrack[i]['mid'])
    Ctracklat.append(CALIOPtrack[i]['lat_mid'])
    if CALIOPtrack[i]['lon_mid']>180:
        Ctracklon.append(CALIOPtrack[i]['lon_mid']-360)
    elif CALIOPtrack[i]['lon_mid']<-180:
        Ctracklon.append(CALIOPtrack[i]['lon_mid']+360)
    else:
        Ctracklon.append(CALIOPtrack[i]['lon_mid'])
        
events=[65,79,93,107,149,178,187,201,215,216,229,230,231,244,259,274,288,317,331,358,411,425,439,453]

# fenetre coupe longitudinale 
ycut=4
widow_lat=10
widow_lon=2
altrangemin=10
altrangemax=26
itera=1
lnwd=5
fs = 16
fig, axs = plt.subplots(1,len(events), sharey=True, gridspec_kw={'wspace': 0.1},figsize=(30,4)) #sharex=True, constrained_layout=True, np.int(len(CALIOPtrack)/ycut)
fig.suptitle('Aerosol scattering ratio 512',fontsize=fs)
axi=0

#for idebut in range(0,len(CALIOPtrack),ycut):
for idebut in events:
    for j in range(len(CALIOPtrack)):
        if CALIOPtrack[j]['cald_index']==idebut:
            break
    ic=j
    ## PLOT CALIOP
    utc = export[ic]['utc']
    alts = export[ic]['alts']
    altii=np.where(alts<=altrangemin)[0][0]
    altiisup=np.where(alts<=altrangemax)[0][0]
    lons = export[ic]['lons']
    lats = export[ic]['lats']
    ilatrangemin=np.where(lats<=Ctracklat[ic]-widow_lat/2)[0][0]
    ilatrangemax=np.where(lats<=Ctracklat[ic]+widow_lat/2)[0][0]
    t512L1_test = export[ic]['t512L1']
    t512L1=t512L1_test[ilatrangemax:ilatrangemin,altiisup:altii]
    t512L1 = np.ma.masked_less(t512L1,0)
    mnd = export[ic]['mnd']
    lbeta512_met = np.log(1000 * mnd * Qs / (kbw*8*np.pi/3))
    malts = export[ic]['malts']
    # calculation of the molecular backscatter
    lbeta512_lid = np.empty(shape=t512L1_test.shape)
    for jy in range(len(lats)):
        lbeta512_lid[jy,:] = np.interp(alts,malts[::-1],lbeta512_met[jy,::-1])
    if medfilter:
        sr512raw = t512L1/np.exp(lbeta512_lid[ilatrangemax:ilatrangemin,altiisup:altii])
        sr512= ss.medfilt(sr512raw,kernel_size=(widthfilter,1))
    else:
        sr512 = t512L1/np.exp(lbeta512_lid[ilatrangemax:ilatrangemin,altiisup:altii])
    im=axs[axi].pcolormesh(lats[ilatrangemax:ilatrangemin],alts[altiisup:altii],sr512.T,cmap=cmap,vmin=vmin,vmax=vmax)

    plt.ylabel('altitude (km)',fontsize=fs)
    fsl=fs-6
    lon_mid=CALIOPtrack[ic]['lon_mid']
    if CALIOPtrack[ic]['date'].month==8:
        if lon_mid>180:
            axs[axi].set_title(str(CALIOPtrack[ic]['date'].day)+'/08 - '+str(np.int(np.abs(lon_mid-360)))+'°W',fontsize=fsl)
        else:
            axs[axi].set_title(str(CALIOPtrack[ic]['date'].day)+'/08 - '+str(np.int(lon_mid))+'°E',fontsize=fsl)
    elif CALIOPtrack[ic]['date'].month==9:
        if lon_mid>180:
            axs[axi].set_title(str(CALIOPtrack[ic]['date'].day)+'/09 - '+str(np.int(np.abs(lon_mid-360)))+'°W',fontsize=fsl)
        else:
            axs[axi].set_title(str(CALIOPtrack[ic]['date'].day)+'/09 - '+str(np.int(lon_mid))+'°E',fontsize=fsl)
    elif CALIOPtrack[ic]['date'].month==10:
        if lon_mid>180:
            axs[axi].set_title(str(CALIOPtrack[ic]['date'].day)+'/10 - '+str(np.int(np.abs(lon_mid-360)))+'°W',fontsize=fsl)
        else:
            axs[axi].set_title(str(CALIOPtrack[ic]['date'].day)+'/10 - '+str(np.int(lon_mid))+'°E',fontsize=fsl)
    plt.xlim(lats[ilatrangemin],lats[ilatrangemax])
    plt.ylim(alts[altii]+0.2,alts[altiisup]-0.2)#alts[altiisup])
    axi+=1
    
for ax in axs:
    ax.label_outer()
    ax.tick_params(labelsize=12)
    ax.set_xlabel('latitude')
#plt.title(txt+" longitude="+str(lon)+'°',fontsize=fs)
plt.colorbar(im)
fig.add_subplot(111, frame_on=False)
plt.tick_params(labelcolor="none", bottom=False, left=False)
#plt.ylabel('baro altitude (km) \n  ',fontsize=fs)
#plt.xlabel('latitude',fontsize=fs)
axs[0].set_ylabel('altitude (km)',fontsize=fs)
#fig.xlabel('latitude',fontsize=fs)
plt.savefig(savfile+'_lon_contour_tot.png',bbox_inches='tight',dpi=300,format='png')
plt.show()