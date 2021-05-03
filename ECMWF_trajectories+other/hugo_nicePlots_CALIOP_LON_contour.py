#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Noice plots.

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
nametrackmodel1 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_A.pkl'
nametrackmodel2 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_B2.pkl'
#nametrackmodel3 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_B2.pkl'
savfile='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/cassure/CONTOUR/'
chemincald=savfile+'selCaliop_'+initials+'.pkl'
nameexport=savfile+'export'
variable='PV'
Qs = 5.167*1.e-31
kbw = 1.0313
#yinf, ysup = 9, 24 
medfilter=True
widthfilter = 21
vmax = 15
vmin=1.5
#lvlcon=[-4:-10]
lvlcon=np.arange(-20,-4.5,15)
cmap = 'gist_ncar'
# load dictionary of traces
with gzip.open(chemincald,'rb') as f:
    box,idate0,nbday,Cald = pickle.load(f)
with gzip.open(nameexport,'rb') as f:
    export=pickle.load(f)
# fenetre coupe longitudinale 
widow_lat=20
widow_lon=2
altrangemin=11
altrangemax=27
itera=1
lnwd=5
fs = 16
track_model1 = pickle.load(open(nametrackmodel1,'rb'))
track_model2= pickle.load(open(nametrackmodel2,'rb'))
#track_model3 = pickle.load(open(nametrackmodel3,'rb'))
events=[30,45,52,67,98] #cassure
#events=[30,45] #cassure
#events=[18,30,97,121,279,319]#30,254

fig, axs = plt.subplots(1,len(events), sharex=True, sharey=True, gridspec_kw={'wspace': 0.1},figsize=(11,4),constrained_layout=True)
fig.suptitle('Contour of anomaly of potential vorticity\n compared to aerosol scattering ratio 512',fontsize=fs)
axi=0

for idebut in events:
    if idebut==67:
        date00=track_model1['dates'][idebut] 
        altitude=track_model1['alts'][idebut]
        masked_alt_min=altrangemin
        if track_model1['lons'][idebut]<0:
            lon_mid=track_model1['lons'][idebut]+360
        else:
            lon_mid=track_model1['lons'][idebut]
        lat_bulle=track_model1['lats'][idebut]
    else:
        date00=track_model2['dates'][idebut] 
        altitude=track_model2['alts'][idebut]
        masked_alt_min=altrangemin
        if track_model2['lons'][idebut]<0:
            lon_mid=track_model2['lons'][idebut]+360
        else:
            lon_mid=track_model2['lons'][idebut]
        lat_bulle=track_model2['lats'][idebut]
    # latrangemin=lat_bulle-widow_lat*0.5
    # latrangemax=lat_bulle+widow_lat*0.5
    lonrangemin=lon_mid-widow_lon/2
    lonrangemax=lon_mid+widow_lon/2
    latrangemin=35
    latrangemax=65
    
    namedats=savfile+'CONTOUR_dats_'+str(idebut)
    if os.path.isfile(namedats):
        with gzip.open(namedats,'rb') as f:
            dats00=pickle.load(f)
    else:    
        dat = ECMWF('FULL-EA',date00,exp='VOZ')
        dat._get_var('O3')
        dat._get_var('T')
        dat._get_var('U') 
        dat._get_var('V')
        dat._get_var('VO')
        dat._get_var('P')
        dat._get_var('PT')
        #dat._get_var('LNSP')
        #dat._get_var('Z')
        #dat._mkpz()
        dat._mkpscale() 
        dat._mkzscale()
        dat._mkp() 
        dat._mkthet() 
        dat._mkpv()
        #dat._mkz() 
        PVmean = np.mean(dat.var['PV'],axis=2)
        dat.var['PVano'] = dat.var['PV'] - PVmean[:,:,np.newaxis]
        VOmean = np.mean(dat.var['VO'],axis=2)
        dat.var['VOano'] = dat.var['VO'] - VOmean[:,:,np.newaxis]
        Pmean=np.mean(dat.var['P'],axis=(1,2))
        #dat.var['Pmean']=Pmean
        #dat.attr['pzscale']=np.log(Pmean/100000)*(-7400)   
        #datp=dat.shift2west(250) 
        dats00 = dat.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))    
        dats00.attr['pzscale']=np.log(Pmean/100000)*(-7.4)   
        ## Enregistrement des données à réutiliser pour ce plot-ci : 
        with gzip.open(namedats,'wb') as f:
           pickle.dump(dats00,f)
    ic=axi
    ## PLOT CALIOP
    utc = export[ic]['utc']
    alts = export[ic]['alts']
    altii=np.where(alts<=masked_alt_min)[0][0]
    altiisup=np.where(alts<=altrangemax)[0][0]
    lons = export[ic]['lons']
    lats = export[ic]['lats']
    ilatrangemin=np.where(lats<=latrangemin)[0][0]
    ilatrangemax=np.where(lats<=latrangemax)[0][0]
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

    #### PLOT CONTOUR
    lon=lon_mid
    # BARO ALT
    #ialtrangemin=np.where(dats00.attr['zscale']<=alts[altii])[0][0] #altrangemin 
    #ialtrangemax=np.where(dats00.attr['zscale']<=alts[altiisup])[0][0] #altrangemax
    # ALTITUDE RÉELLE
    ialtrangemin=np.where(dats00.attr['pzscale']<=alts[altii])[0][0] #altrangemin
    ialtrangemax=np.where(dats00.attr['pzscale']<=alts[altiisup])[0][0] #altrangemax
    levs=(ialtrangemax,ialtrangemin)
    scale=10**6
    try:
        pos=np.where(dats00.attr['lons']>=lon)[0][0]
        print('pos',pos)
    except:
        print('lon out of range')
    if levs[0]==None: l1=29
    else: l1 = levs[0]
    if levs[1]==None: l2=115
    else: l2 = levs[1]
    #if Z variable is available, lets use it, if not use zscale
    try:
        zz1 = 0.5*(dats00.var['Z'][l1-1:l2+1,pos, :] + dats00.var['Z'][l1:l2+2,pos,:])/1000
        zz = np.empty((zz1.shape[0],zz1.shape[1]+1))
        zz[:,1:-1] = 0.5*(zz1[:,1:]+zz1[:,:-1])
        zz[:,0] = zz1[:,0]
        zz[:,-1] = zz1[:,-1]
        print(zz.shape,len(lons))
        #axs[axi].pcolormesh(lons,zz,scale*dats00.var[var][l1:l2+1,pos,:],vmin=clim[0],vmax=clim[1],cmap=cmap)
        axs[axi].contour(scale*dats00.var[variable+'ano'][l1:l2+1,:,pos],levels=lvlcon,extent=(dats00.attr['lats'][0], dats00.attr['lats'][-1],dats00.attr['Z'][ialtrangemax],dats00.attr['Z'][ialtrangemin]-2),colors='red',linewidths=2.5)#extent=(lats[ilatrangemin],lats[ilatrangemax],alts[altiisup],alts[altii])#latrangemin,latrangemax,altrangemin,altrangemax
        plt.ylabel('altitude (km)',fontsize=fs)
        print('USE Z')
    except:
        axs[axi].contour(scale*dats00.var[variable+'ano'][l1:l2+1,:,pos],levels=lvlcon,extent=(dats00.attr['lats'][0], dats00.attr['lats'][-1],dats00.attr['pzscale'][ialtrangemax]-2.2,dats00.attr['pzscale'][ialtrangemin]-2.2),colors='red',linewidths=2.5)#extent=(lats[ilatrangemin],lats[ilatrangemax],alts[altiisup],alts[altii])#latrangemin,latrangemax,altrangemin,altrangemax

        plt.ylabel('altitude (km)',fontsize=fs)

    #axs[axi].scatter(lat_bulle,altitude,linewidth=3,color='red',marker='+',s=20**2)
    axs[axi].scatter(lat_bulle,dats00.attr['pzscale'][np.where(dats00.attr['zscale']<=altitude)[0][0]]-1,linewidth=3,color='red',marker='+',s=20**2)
    if track_model2['dates'][idebut].month==8:
        if lon_mid>180:
            axs[axi].set_title(str(track_model2['dates'][idebut].day)+' Aug    '+str(np.int(np.abs(lon_mid-360)))+'°W',fontsize=fs-3)
        else:
            axs[axi].set_title(str(track_model2['dates'][idebut].day)+' Aug    '+str(np.int(lon_mid))+'°E',fontsize=fs-3)
    elif track_model2['dates'][idebut].month==9:
        if lon_mid>180:
            axs[axi].set_title(str(track_model2['dates'][idebut].day)+' Sept    '+str(np.int(np.abs(lon_mid-360)))+'°W',fontsize=fs-3)
        else:
            axs[axi].set_title(str(track_model2['dates'][idebut].day)+' Sept    '+str(np.int(lon_mid))+'°E',fontsize=fs-3)
    plt.xlim(lats[ilatrangemin],lats[ilatrangemax])
    plt.ylim(alts[altii]+0.2,21.5)#alts[altiisup])
    axi+=1
    
#for ix in range(len(axs)-1):
    #axs[ix].plot([lats[ilatrangemax]-0.1,lats[ilatrangemax]-0.1],[alts[altiisup],alts[altii]],linewidth=9,color='white')
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
plt.savefig(savfile+str(variable)+'_lon_contour_tot.png',bbox_inches='tight',dpi=300,format='png')
plt.show()