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
variable='O3'
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
else:    
    dat = ECMWF('FULL-EA',date00,exp='VOZ') 
    dat._get_var('O3')
    dat._get_var('T')
    dat._get_var('U') 
    dat._get_var('V')
    dat._get_var('VO')
    dat._get_var('P')
    dat._get_var('PT')
    dat._mkpscale() 
    dat._mkzscale()
    dat._mkp() 
    dat._mkthet() 
    dat._mkpv()
    PVmean = np.mean(dat.var[variable],axis=2)
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    #datp=dat.shift2west(260) 
    dats00 = dat.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
        
    ## Enregistrement des données à réutiliser pour ce plot-ci : 
    with gzip.open(savfile+'dats00_'+str(idebut),'wb') as f:
        pickle.dump(dats00,f)

        
## map de données :
ialts=np.where(dats00.attr['zscale']<=altitude)[0][0] 
jy = np.where(dats00.attr['lats']>=lat_bulle)[0][0]
events=[170,195,218,233,254]
#events=np.arange(98,149,10)
#jourdelta=4
#itrack=0
#date00 = date00 + timedelta(days=jourdelta)
#while  date00.day!=25:
#while itrack!=254:
    # for itrack in range(idebut,255):
    #     if track_model3['dates'][itrack].day==date00.day and track_model3['dates'][itrack].hour==date00.hour:
    #         break
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
    else:    
        dat = ECMWF('FULL-EA',date00,exp='VOZ') 
        dat._get_var('O3')
        dat._get_var('T')
        dat._get_var('U') 
        dat._get_var('V')
        dat._get_var('VO')
        dat._get_var('P')
        dat._get_var('PT')
        dat._mkpscale() 
        dat._mkzscale()
        dat._mkp() 
        dat._mkthet() 
        dat._mkpv()
        PVmean = np.mean(dat.var[variable],axis=2)
        dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
        dats = dat.extract(varss='All',lonRange=(lon_min_w,lon_max_w),latRange=(lat_min_w,lat_max_w))
        with gzip.open(namedate00,'wb') as f:
            pickle.dump(dats,f)
    # PVmean = np.mean(dats.var[variable],axis=2)
    # dats.var[variable+'ano'] = dats.var[variable] - PVmean[:,:,np.newaxis]
    kz1 = np.where(dats00.attr['zscale']<=track_model3['alts'][itrack])[0][0]
    jy2 = np.where(dats.attr['lats']==track_model3['lats'][itrack])[0][0]
    dats00.var[variable+'ano'][ialts,jy1-int(round(fenetre_lat/2)):jy1+int(round(fenetre_lat/2)),ix1-int(round(fenetre_lon/2)):ix1+int(round(fenetre_lon/2))+1] = dats.var[variable+'ano'][kz1,:,:]
    dats00.var[variable+'anoZ'][alt_max_w:alt_min_w,jy,ix1-int(round(fenetre_lon/2)):ix1+int(round(fenetre_lon/2))+1]=dats.var[variable+'ano'][alt_max_w:alt_min_w,jy2,:]
    #date00 = date00 + timedelta(days=jourdelta)

#dats00.show(variable+'ano',ialts,cLines=False,projec=None,txt=str(variable)+' alt')   

#with gzip.open(savfile+'dats00_tot'),'wb') as f:
    #pickle.dump(dats00,f)

jet3= plt.get_cmap('binary')
colors3 = iter(jet3(np.linspace(0,1,colori))) 
colors32 = iter(jet3(np.linspace(0,1,colori))) 
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
lnwd=3

#### PLOT ZONAL
lat=lat_bulle
ialtrangemin=np.where(dats.attr['zscale']<=altrangemin)[0][0]
ialtrangemax=np.where(dats.attr['zscale']<=altrangemax)[0][0]
levs=(ialtrangemax,ialtrangemin)
log=False
#clim=(None,None)
cmap=mymap
if variable=='PV':
    txt='Anomaly of potential vorticity ($10^{-6}$ m$^{2}$ s$^{-1}$ K kg$^{-1}$)'
    clim=(-15,15)
    scale=10**6
elif variable=='O3':
    txt='Anomaly of ozone ($10^{-6}$ kg/kg)'
    clim=(0,-1)
    scale=10**6
try:
    pos=np.where(dats00.attr['lats']>=lat)[0][0]
    print('pos',pos)
except:
    print('lat out of range')
fig2 = plt.figure(figsize=(11,4))
fs = 15
ax2 = fig2.add_subplot(111)
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
    iax2=ax2.pcolormesh(lons,zz,scale*dats00.var[variable+'anoZ'][l1:l2+1,pos,:],
                    vmin=clim[0],vmax=clim[1],cmap=cmap)
    plt.ylabel('altitude (km)',fontsize=fs)
    print('USE Z')
except(KeyError):
    iax2=ax2.pcolormesh(lons,dats00.attr['zscale_i'][l1:l2+2],scale*dats00.var[variable+'anoZ'][l1:l2+1,pos, :],
                vmin=clim[0],vmax=clim[1],cmap=cmap)

    plt.ylabel('baro altitude (km)',fontsize=fs)
ax2.tick_params(labelsize=16)
plt.xlabel('longitude',fontsize=fs)

plt.title(txt+" latitude="+str(lat)+'°',fontsize=fs)
cbar = fig2.colorbar(iax2)
cbar.ax.tick_params(labelsize=fs)

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
    ax2.plot([lon_min_w,lon_min_w],[alt_min_w,alt_max_w],linewidth=2,color='white') # gauche
    ax2.plot([lon_max_w,lon_max_w],[alt_min_w,alt_max_w],linewidth=2,color='white') # droite
    ax2.plot([lon_min_w,lon_max_w],[alt_max_w,alt_max_w],linewidth=2,color='white') # haut
    ax2.plot([lon_min_w,lon_max_w],[alt_min_w,alt_min_w],linewidth=2,color='white') # bas
    xleg2,yleg2=proj2legend(track_model3['lons'][itrack]-9.2,track_model3['alts'][itrack]-altrangemin-1.9,lonrangemin,lonrangemax,altrangemin,altrangemax)
    #ax.legend(loc=(xleg,yleg))
    if track_model3['dates'][itrack].month==8:
        textstr2=str(track_model3['dates'][itrack].day)+' Aug'
    elif track_model3['dates'][itrack].month==9:
        textstr2=str(track_model3['dates'][itrack].day)+' Sept'
    ax2.text(xleg2,yleg2,textstr2, transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=props2)

for i in range(idebut,263,itera):
    if  track_model3['dates'][i].month==9 and 3<track_model3['dates'][i].day<16:
        ax2.plot([lons3[i-itera],lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color='white',marker='x',alpha=0.5)
    else:
        ax2.plot([lons3[i-itera],lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color='white',alpha=0.5)

plt.xticks([30,60,90,120,150,180],['30°E','60°E','90°E','120°E','150°E','180°E'])
plt.savefig(savfile+str(variable)+'_lat.png',bbox_inches='tight',dpi=300,format='png')
plt.show()