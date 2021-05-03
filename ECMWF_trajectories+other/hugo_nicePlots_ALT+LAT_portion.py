#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Noice plots.

@author: Hugo Lestrelin

"""
import sys,os  
import socket      
if 'ciclad' in socket.gethostname():
    import sys,os
    sys.path.append('/home/hlestrel/pylib')                
from datetime import datetime,timedelta
#from ECMWF_N_test_proj import ECMWF 
from ECMWF_N import ECMWF 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm
import pickle,gzip
from matplotlib import gridspec
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from cartopy import feature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
from progress_bar import printProgressBar
from mpl_toolkits.mplot3d import Axes3D
import gc

listcolors=['#161d58','#253494','#2850a6','#2c7fb8','#379abe','#41b6c4',
            '#71c8bc','#a1dab4','#d0ecc0','#ffffcc','#fef0d9','#fedeb1',
            '#fdcc8a','#fdac72','#fc8d59','#ef6b41','#e34a33','#cb251a',
            '#b30000','#7f0000']
# homemade color table
listcolors_sw=[listcolors[1],listcolors[0],listcolors[3],listcolors[2],\
               listcolors[5],listcolors[4],listcolors[7],listcolors[6],\
               listcolors[9],listcolors[8],listcolors[11],listcolors[10],\
               listcolors[13],listcolors[12],listcolors[15],listcolors[14],\
               listcolors[17],listcolors[16],listcolors[19],listcolors[18]]
mymap=colors.ListedColormap(listcolors)
mymap_sw=colors.ListedColormap(listcolors_sw)
#mymap='RdGy'

## fonction transform proj -> legend
def proj2legend(x,y,x0,x1,y0,y1):
    xproj=x1-x0
    yproj=y1-y0
    xlegend=x/xproj
    ylegend=y/yproj
    return(xlegend,ylegend)

initials='CAN2017'
nametrackmodel3 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_B2.pkl'
savfile='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/'
variable='PV'
fenetre_lon=19 # fenêtre des subplots
fenetre_lat=15
fenetre_alt=4
# fenetre image fond
lonrangemin=0
lonrangemax=200
latrangemin=17
latrangemax=85
altrangemin=14
altrangemax=23.5


projplate = ccrs.PlateCarree(central_longitude=180)
proj=projplate 
figsize=(20,6)
fs=15
fig = plt.figure(figsize=figsize)#,constrained_layout=True)
fontsize=fs


###### UPDATE TRACK
colori=197
idebut=71
track_model3 = pickle.load(open(nametrackmodel3,'rb'))
idx3=len(track_model3['dates'])
date00=track_model3['dates'][idebut] #choix de l'indice pour l'image de fond
altitude=track_model3['alts'][idebut]
lon_mid=track_model3['lons'][idebut]
lat_bulle=track_model3['lats'][idebut]
if (lon_mid+fenetre_lon/2)<359:
    loninf = lon_mid-fenetre_lon/2
    lonsup = lon_mid+fenetre_lon/2
elif (lon_mid+fenetre_lon/2)>359:
    loninf = lon_mid-fenetre_lon/2-359
    lonsup = lon_mid+fenetre_lon/2-359
if (lat_bulle-fenetre_lat/2)>-85 and (lat_bulle+fenetre_lat/2)<85:
    latinf = lat_bulle-fenetre_lat/2
    latsup = lat_bulle+fenetre_lat/2
elif (lat_bulle-fenetre_lat/2)<-85:
    latinf = -85
    latsup = lat_bulle+fenetre_lat/2
elif (lat_bulle+fenetre_lat/2)>85:
    latinf = lat_bulle-fenetre_lat/2
    latsup = 85
    
if os.path.isfile(savfile+'dats00_'+str(idebut)):
    with gzip.open(savfile+'dats00_'+str(idebut),'rb') as f:
        dats00=pickle.load(f)
        
## map de données :
ialts=np.where(dats00.attr['zscale']<=altitude)[0][0] 
jy = np.where(dats00.attr['lats']>=lat_bulle)[0][0]
events=[170,195,218,233,254]
## double des données pour la coupe zonale :
PVmean = np.mean(dats00.var[variable],axis=2)
dats00.var[variable+'ano'] = dats00.var[variable] - PVmean[:,:,np.newaxis]
dats00.var[variable+'anoZ'] = dats00.var[variable+'ano'].copy()
for itrack in events:
    date00=track_model3['dates'][itrack]
    namedate00=savfile+'dats00_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
    lon_min_w=track_model3['lons'][itrack]-fenetre_lon/2
    lon_max_w=track_model3['lons'][itrack]+fenetre_lon/2
    lat_min_w=track_model3['lats'][itrack]-fenetre_lat/2
    lat_max_w=track_model3['lats'][itrack]+fenetre_lat/2
    alt_min_w=np.where(dats00.attr['zscale']<=track_model3['alts'][itrack]-fenetre_alt*0.4)[0][0]
    alt_max_w=np.where(dats00.attr['zscale']<=track_model3['alts'][itrack]+fenetre_alt*0.6)[0][0]
    ix1 = int(track_model3['lons'][itrack]-dats00.attr['lons'][0])
    jy1 = int(track_model3['lats'][itrack]-dats00.attr['lats'][0])
    if os.path.isfile(namedate00):
        with gzip.open(namedate00,'rb') as f:
            dats=pickle.load(f)
    #PVmean = np.mean(dats.var[variable],axis=2)
    #dats.var[variable+'ano'] = dats.var[variable] - PVmean[:,:,np.newaxis]
    kz1 = np.where(dats00.attr['zscale']<=track_model3['alts'][itrack])[0][0]
    jy2 = np.where(dats.attr['lats']==track_model3['lats'][itrack])[0][0]
    dats00.var[variable+'ano'][ialts,jy1-int(round(fenetre_lat/2)):jy1+int(round(fenetre_lat/2)),ix1-int(round(fenetre_lon/2)):ix1+int(round(fenetre_lon/2))+1] = dats.var[variable+'ano'][kz1,:,:]
    dats00.var[variable+'anoZ'][alt_max_w:alt_min_w,jy,ix1-int(round(fenetre_lon/2)):ix1+int(round(fenetre_lon/2))+1]=dats.var[variable+'ano'][alt_max_w:alt_min_w,jy2,:]
 
lons3=[]
lats3=[]
alts3=[]
for i in range(idx3):
    alts3.append(track_model3['alts'][i])
    lats3.append(track_model3['lats'][i])
    if track_model3['lons'][i]>259:
        lons3.append(track_model3['lons'][i]-360)
    elif track_model3['lons'][i]<-180:
        lons3.append(track_model3['lons'][i]+360)
    else:
        lons3.append(track_model3['lons'][i])
itera=1
lnwd=2

ax0 = plt.subplot(211,projection = proj)
##### PLOT horizontal
lev=ialts
cardinal_level=True
clim=(-15,15)
cmap=mymap
scale=10**6
aspect=1
if len(dats00.var[variable+'ano'].shape) == 3:
    if (cardinal_level==False) | (lev > dats00.nlev-1):
        clev = np.abs(np.array(dats00.attr['levs'])-lev).argmin()
    else:
        clev = lev
    buf = dats00.var[variable+'ano'][clev,:,:]
else:
    buf = dats00.var[variable+'ano']
fs=15
cm_lon =0
ctrl_lon=(dats00.attr['lons'][0]+dats00.attr['lons'][-1])/2
ctrl_lat=(dats00.attr['lats'][0]+dats00.attr['lats'][-1])/2
if dats00.attr['lons'][-1] > 180: cm_lon=180
iax = ax0.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats00.attr['lons'][0]-cm_lon, dats00.attr['lons'][-1]-cm_lon,dats00.attr['lats'][0], dats00.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
xlocs = None
if cm_lon == 180:
        interx = 30
        minx = dats00.attr['lons'][0] + interx - dats00.attr['lons'][0]%interx
        xlocs = list(np.arange(minx,181,interx))+list(np.arange(interx-180,dats00.attr['lons'][-1]-360,interx))
gl = ax0.gridlines(draw_labels=True, xlocs=xlocs,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax0.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax0.coastlines('50m')
gl.top_labels = False
gl.bottom_labels = False
gl.right_labels = False
gl.left_labels = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
#gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs} 
props = dict(boxstyle='round', facecolor='white', alpha=0.8)
## BOX HORIZ
for itrack in events:
    lon_min_w=track_model3['lons'][itrack]-fenetre_lon/2-1
    lon_max_w=track_model3['lons'][itrack]+fenetre_lon/2+1
    lat_min_w=track_model3['lats'][itrack]-fenetre_lat/2-1
    lat_max_w=track_model3['lats'][itrack]+fenetre_lat/2
    ax0.plot([lon_min_w-180,lon_min_w-180],[lat_min_w,lat_max_w],linewidth=2,color='black') # gauche
    ax0.plot([lon_max_w-180,lon_max_w-180],[lat_min_w,lat_max_w],linewidth=2,color='black') # droite
    ax0.plot([lon_min_w-180,lon_max_w-180],[lat_max_w,lat_max_w],linewidth=2,color='black') # haut
    ax0.plot([lon_min_w-180,lon_max_w-180],[lat_min_w,lat_min_w],linewidth=2,color='black') # bas
    xleg,yleg=proj2legend(track_model3['lons'][itrack]-8,track_model3['lats'][itrack]-latrangemin-11,lonrangemin,lonrangemax,latrangemin,latrangemax)
    #ax.legend(loc=(xleg,yleg))
    if track_model3['dates'][itrack].month==8:
        textstr=str(track_model3['dates'][itrack].day)+' Aug'
    elif track_model3['dates'][itrack].month==9:
        textstr=str(track_model3['dates'][itrack].day)+' Sept'
    ax0.text(xleg,yleg,textstr, transform=ax0.transAxes, fontsize=14, verticalalignment='top', bbox=props)

## TRACK
for i in range(idebut,263,itera):
    if  track_model3['dates'][i].month==9 and 3<track_model3['dates'][i].day<16:
        ax0.plot([lons3[i-itera]-180,lons3[i]-180],[lats3[i-itera],lats3[i]],linewidth=lnwd,color='k',marker='x',alpha=0.5)
    else:
        ax0.plot([lons3[i-itera]-180,lons3[i]-180],[lats3[i-itera],lats3[i]],linewidth=lnwd,color='k')
#plt.xlim(0,200)

#### PLOT ZONAL
ax1 = plt.subplot(212,sharex=ax0)
lat=lat_bulle
ialtrangemin=np.where(dats.attr['zscale']<=altrangemin)[0][0]
ialtrangemax=np.where(dats.attr['zscale']<=altrangemax)[0][0]
levs=(ialtrangemax,ialtrangemin)
txt='Anomaly of potential vorticity ($10^{-6}$ m$^{2}$ s$^{-1}$ K kg$^{-1}$)'
log=False
clim=(-15,15)
#clim=(None,None)
cmap=mymap
scale=10**6
try:
    pos=np.where(dats00.attr['lats']>=lat)[0][0]
    print('pos',pos)
except:
    print('lat out of range')
fs = 15
if levs[0]==None: l1=29
else: l1 = levs[0]
if levs[1]==None: l2=115
else: l2 = levs[1]
lons = np.arange(dats00.attr['lons'][0]-0.5*dats00.attr['dlo'],dats00.attr['lons'][-1]+dats00.attr['dlo'],dats00.attr['dlo'])
try:
    zz1 = 0.5*(dats00.var['Z'][l1-1:l2+1,pos, :] + dats00.var['Z'][l1:l2+2,pos,:])/1000
    zz = np.empty((zz1.shape[0],zz1.shape[1]+1))
    zz[:,1:-1] = 0.5*(zz1[:,1:]+zz1[:,:-1])
    zz[:,0] = zz1[:,0]
    zz[:,-1] = zz1[:,-1]
    print(zz.shape,len(lons))
    iax2=ax1.pcolormesh(lons,zz,scale*dats00.var[variable+'anoZ'][l1:l2+1,pos,:],
                    vmin=clim[0],vmax=clim[1],cmap=cmap)
    plt.ylabel('altitude (km)',fontsize=fs)
    print('USE Z')
except(KeyError):
    iax2=ax1.pcolormesh(lons,dats00.attr['zscale_i'][l1:l2+2],scale*dats00.var[variable+'anoZ'][l1:l2+1,pos, :],
                vmin=clim[0],vmax=clim[1],cmap=cmap)

    plt.ylabel('baro altitude (km)',fontsize=fs)
ax1.tick_params(labelsize=16)
plt.xlabel('longitude',fontsize=fs)
ax1.set_xlim(0,200)

props2 = dict(boxstyle='round', facecolor='white', alpha=0.8)
## BOX ZONAL
for itrack in events:
    lon_min_w=track_model3['lons'][itrack]-fenetre_lon/2-1
    lon_min_wi=np.where(dats00.attr['lons']>=lon_min_w)[0][0]
    lon_max_w=track_model3['lons'][itrack]+fenetre_lon/2+1
    lon_max_wi=np.where(dats00.attr['lons']>=lon_max_w)[0][0]
    alt_min_w=track_model3['alts'][itrack]-fenetre_alt*0.4
    alt_max_w=track_model3['alts'][itrack]+fenetre_alt*0.6-0.1
    #ax2.contour(scale*dats00.var[variable+'anoZ'][l1:l2+1,pos,lon_min_wi:lon_max_wi],extent=(lon_min_w, lon_max_w,alt_min_w,alt_max_w),levels=[-10],origin='lower')
    #ax2.contour(scale*dats00.var[variable+'anoZ'][l1:l2+1,pos,:],extent=(lonrangemin, lonrangemax,altrangemax,altrangemin-1.8),levels=[-10])#,origin='lower')
    ax1.plot([lon_min_w,lon_min_w],[alt_min_w,alt_max_w],linewidth=2,color='black') # gauche
    ax1.plot([lon_max_w,lon_max_w],[alt_min_w,alt_max_w],linewidth=2,color='black') # droite
    ax1.plot([lon_min_w,lon_max_w],[alt_max_w,alt_max_w],linewidth=2,color='black') # haut
    ax1.plot([lon_min_w,lon_max_w],[alt_min_w,alt_min_w],linewidth=2,color='black') # bas
    xleg2,yleg2=proj2legend(track_model3['lons'][itrack]-9.2,track_model3['alts'][itrack]-altrangemin-1.9,lonrangemin,lonrangemax,altrangemin,altrangemax)
    #ax.legend(loc=(xleg,yleg))
    if track_model3['dates'][itrack].month==8:
        textstr2=str(track_model3['dates'][itrack].day)+' Aug'
    elif track_model3['dates'][itrack].month==9:
        textstr2=str(track_model3['dates'][itrack].day)+' Sept'
    ax1.text(xleg2,yleg2,textstr2, transform=ax1.transAxes, fontsize=12, verticalalignment='top', bbox=props2)

for i in range(idebut,263,itera):
    if  track_model3['dates'][i].month==9 and 3<track_model3['dates'][i].day<16:
        ax1.plot([lons3[i-itera],lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color='k',marker='x',alpha=0.5)
    else:
        ax1.plot([lons3[i-itera],lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color='k')

fig.subplots_adjust(hspace=0.1)#,wspace=0.5,top=0.925,left=0.)
pos_cax = fig.add_axes([1.05, 0.15, 0.05, 0.3])
cbar=fig.colorbar(iax,cax=pos_cax)
cbar.ax.tick_params(labelsize=fs) 
fig.suptitle('Anomaly of potential vorticity ($10^{-6}$ m$^{2}$ s$^{-1}$ K kg$^{-1}$)',fontsize=fs)
plt.savefig(savfile+str(variable)+'_LAT+ALT.png',bbox_inches='tight',dpi=300,format='png')
plt.show()