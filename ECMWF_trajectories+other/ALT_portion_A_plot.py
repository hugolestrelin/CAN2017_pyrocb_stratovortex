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
nametrackmodel3 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_A.pkl'
savfile='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/VA_portion/'
variable='O3'
fenetre_lon=19 # fenêtre des subplots
fenetre_lat=15
fenetre_alt=4
# fenetre image fond
lonrangemin=-180
lonrangemax=100
latrangemin=15
latrangemax=75
altrangemin=14
altrangemax=23.5
###### UPDATE TRACK
idebut=12
[trac_PV,track_model3] = pickle.load(open(savfile+'VortexA-track.pkl','rb'))
idx3=len(track_model3['dates'])
colori=idx3
date00=track_model3['dates'][idebut] #choix de l'indice pour l'image de fond
altitude=track_model3['z'][idebut]
lon_mid=track_model3['lons'][idebut]
lat_bulle=track_model3['lats'][idebut]
    
namedats=savfile+'dats00_'+str(idebut)
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
    dat._mkpscale() 
    dat._mkzscale()
    dat._mkp() 
    dat._mkthet() 
    dat._mkpv()
    PVmean = np.mean(dat.var['PV'],axis=2)
    dat.var['PVano'] = dat.var['PV'] - PVmean[:,:,np.newaxis]
    O3mean = np.mean(dat.var['O3'],axis=2)
    dat.var['O3ano'] = dat.var['O3'] - O3mean[:,:,np.newaxis]
    Tmean = np.mean(dat.var['T'],axis=2)
    dat.var['Tano'] = dat.var['T'] - Tmean[:,:,np.newaxis]
    datp=dat.shift2west(140) 
    dats00 = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
        
    ## Enregistrement des données à réutiliser pour ce plot-ci : 
    with gzip.open(namedats,'wb') as f:
        pickle.dump(dats00,f)

if os.path.isfile(savfile+'dats_meantotO3'):
    with gzip.open(namedats,'rb') as f:
        dats00=pickle.load(f)
else:
    for i in range(1,222):
        namedats1=savfile+'min_O3/dats_'+str(i)
        print(namedats1)
        with gzip.open(namedats1,'rb') as f:
            datsi=pickle.load(f)
        dats00.var['O3ano'][:,:,0:281]+=datsi.var['O3ano']
    dats00.var['O3ano']=dats00.var['O3ano']/222
    with gzip.open(savfile+'dats_meantotO3','wb') as f:
        pickle.dump(dats00,f)
        
## map de données :
ialts=np.where(dats00.attr['zscale']<=altitude)[0][0] 
jy = np.where(dats00.attr['lats']>=lat_bulle)[0][0]
#dats00.var[variable+'ano'][ialts,:,:]='NaN' #sans image de fond
#events=[98,170,195,218,233,254]
events=[12,28,64,96,112,136,160,176,192] # [12,28,68,92,108,132,156,172,192]
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
for itrack in events:
    date00=track_model3['dates'][itrack]
    namedate00=savfile+'dats00_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
    lon_min_w=track_model3['lons'][itrack]-fenetre_lon/2
    lon_max_w=track_model3['lons'][itrack]+fenetre_lon/2
    lat_min_w=track_model3['lats'][itrack]-fenetre_lat/2
    lat_max_w=track_model3['lats'][itrack]+fenetre_lat/2
    alt_min_w=np.where(dats00.attr['zscale']<=track_model3['z'][itrack]-fenetre_alt*0.4)[0][0]
    alt_max_w=np.where(dats00.attr['zscale']<=track_model3['z'][itrack]+fenetre_alt*0.6)[0][0]
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
        PVmean = np.mean(dat.var['PV'],axis=2)
        dat.var['PVano'] = dat.var['PV'] - PVmean[:,:,np.newaxis]
        O3mean = np.mean(dat.var['O3'],axis=2)
        dat.var['O3ano'] = dat.var['O3'] - O3mean[:,:,np.newaxis]
        Tmean = np.mean(dat.var['T'],axis=2)
        dat.var['Tano'] = dat.var['T'] - Tmean[:,:,np.newaxis]
        datp=dat.shift2west(140) 
        dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
        with gzip.open(namedate00,'wb') as f:
            pickle.dump(dats,f)
    kz1 = np.where(dats00.attr['zscale']<=track_model3['z'][itrack])[0][0]
    jy2 = np.where(dats.attr['lats']==track_model3['lats'][itrack])[0][0]
    dats00.var[variable+'ano'][ialts,jy1-int(round(fenetre_lat/2)):jy1+int(round(fenetre_lat/2)),ix1-int(round(fenetre_lon/2)):ix1+int(round(fenetre_lon/2))+1] = dats.var[variable+'ano'][kz1,jy1-int(round(fenetre_lat/2)):jy1+int(round(fenetre_lat/2)),ix1-int(round(fenetre_lon/2)):ix1+int(round(fenetre_lon/2))+1]
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
    alts3.append(track_model3['z'][i])
    lats3.append(track_model3['lats'][i])
    if track_model3['lons'][i]>lonrangemax:
        lons3.append(track_model3['lons'][i]-360)
    #elif track_model3['lons'][i]<-180:
        #lons3.append(track_model3['lons'][i]+360)
    else:
        lons3.append(track_model3['lons'][i])
itera=1
lnwd=3
# 
# ##### PLOT horizontal
lev=ialts
cardinal_level=True
if variable=='PV':
    txt='Anomaly of potential vorticity ($10^{-6}$ m$^{2}$ s$^{-1}$ K kg$^{-1}$)'
    clim=(-15,15)
    scale=10**6
elif variable=='O3':
    txt='Anomaly of ozone ($10^{-6}$ kg/kg)'
    clim=(0.75,-0.75)
    scale=10**6
figsize=(11,4)
cmap=mymap
cLines=None
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
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
proj=projplate
fig1 = plt.figure(figsize=figsize)#,constrained_layout=True)
fig1.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
ax = plt.axes(projection = proj)
iax = ax.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats00.attr['lons'][0]-cm_lon, dats00.attr['lons'][-1]-cm_lon,dats00.attr['lats'][0], dats00.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
xlocs = None
if cm_lon == 180:
        interx = 30
        minx = dats00.attr['lons'][0] + interx - dats00.attr['lons'][0]%interx
        xlocs = list(np.arange(minx,181,interx))+list(np.arange(interx-180,dats00.attr['lons'][-1]-360,interx))
gl = ax.gridlines(crs=projplate,draw_labels=True, xlocs=xlocs,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
if cLines is not None:
            ax.contour(dats00.var[variable+'ano'][lev,:,:],transform=projplate,extent=(dats00.attr['lons'][0]-cm_lon,
                        dats00.attr['lons'][-1]-cm_lon,dats00.attr['lats'][0],dats00.attr['lats'][-1]),levels=cLines,origin='lower')
ax.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax.coastlines('50m')
gl.top_labels = False
gl.right_labels = False
gl.left_labels = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
if txt is None:
    plt.title(var+' '+str(lev),fontsize=fs)
else:
    plt.title(txt,fontsize=fs)
axpos = ax.get_position()
pos_x = axpos.x0 + axpos.x0 + axpos.width + 0.01
pos_cax = fig1.add_axes([pos_x,axpos.y0,0.04,axpos.height])
cbar=fig1.colorbar(iax,cax=pos_cax)
cbar.ax.tick_params(labelsize=fs)

#ax=dats00.show('PVano',ialts,projec=None,show=False,txt='Anomaly of potential vorticity ($10^{-6}$ m$^{2}$ s$^{-1}$ K kg$^{-1}$)',clim=(-15,15),scale=10**6,figsize=(11,4),aspect=1)


#track_model3['lons'][itrack]=track_model3['lons'][itrack]+180
props = dict(boxstyle='round', facecolor='white', alpha=0.85)
## BOX HORIZ
for itrack in events:
    lon_min_w=track_model3['lons'][itrack]-fenetre_lon/2-1
    lon_max_w=track_model3['lons'][itrack]+fenetre_lon/2+1
    lat_min_w=track_model3['lats'][itrack]-fenetre_lat/2-0.6
    lat_max_w=track_model3['lats'][itrack]+fenetre_lat/2
    if variable=='O3':
        ax.plot([lon_min_w,lon_min_w],[lat_min_w,lat_max_w],linewidth=2,color='black') # gauche
        ax.plot([lon_max_w,lon_max_w],[lat_min_w,lat_max_w],linewidth=2,color='black') # droite
        ax.plot([lon_min_w,lon_max_w],[lat_max_w,lat_max_w],linewidth=2,color='black') # haut
        ax.plot([lon_min_w,lon_max_w],[lat_min_w,lat_min_w],linewidth=2,color='black') # bas
    elif variable=='PV':
        ax.plot([lon_min_w,lon_min_w],[lat_min_w,lat_max_w],linewidth=2,color='black') # gauche
        ax.plot([lon_max_w,lon_max_w],[lat_min_w,lat_max_w],linewidth=2,color='black') # droite
        ax.plot([lon_min_w,lon_max_w],[lat_max_w,lat_max_w],linewidth=2,color='black') # haut
        ax.plot([lon_min_w,lon_max_w],[lat_min_w,lat_min_w],linewidth=2,color='black') # bas
    xleg,yleg=proj2legend(track_model3['lons'][itrack]+180-6.2,track_model3['lats'][itrack]-latrangemin+11,lonrangemin,lonrangemax,latrangemin,latrangemax)
    #ax.legend(loc=(xleg,yleg))
    if track_model3['dates'][itrack].month<10:
        if track_model3['dates'][itrack].day<10:
            textstr='0'+str(track_model3['dates'][itrack].day)+'/0'+str(track_model3['dates'][itrack].month)#+'/'+str(track_model3['dates'][itrack].year)
        else:
            textstr=str(track_model3['dates'][itrack].day)+'/0'+str(track_model3['dates'][itrack].month)#+'/'+str(track_model3['dates'][itrack].year)
    else:
        if track_model3['dates'][itrack].day<10:
            textstr='0'+str(track_model3['dates'][itrack].day)+'/'+str(track_model3['dates'][itrack].month)#+'/'+str(track_model3['dates'][itrack].year)
        else:
            textstr=str(track_model3['dates'][itrack].day)+'/'+str(track_model3['dates'][itrack].month)#+'/'+str(track_model3['dates'][itrack].year)
    ax.text(xleg,yleg,textstr, transform=ax.transAxes, fontsize=9, verticalalignment='top', bbox=props)

## TRACK
#for i in range(idebut,263,itera):
if variable=='O3':
    for i in range(1,len(track_model3['dates'])-1):
        if  track_model3['dates'][i].month==9 and 3<track_model3['dates'][i].day<16:
            ax.plot([lons3[i-itera],lons3[i]],[lats3[i-itera],lats3[i]],linewidth=3,color='navy')#,marker='x',alpha=1)
        else:
            ax.plot([lons3[i-itera],lons3[i]],[lats3[i-itera],lats3[i]],linewidth=3,color='navy',alpha=1)
elif variable=='PV':
    for i in range(1,len(track_model3['dates'])-1):
        if  track_model3['dates'][i].month==9 and 3<track_model3['dates'][i].day<16:
            ax.plot([lons3[i-itera],lons3[i]],[lats3[i-itera],lats3[i]],linewidth=1,color='navy')#,marker='x',alpha=1)
        else:
            ax.plot([lons3[i-itera],lons3[i]],[lats3[i-itera],lats3[i]],linewidth=1,color='navy',alpha=1)

#fig.add_subplot(111, frame_on=False)
#plt.tick_params(labelcolor="none", bottom=False, left=False)
#plt.ylabel('baro altitude (km) \n  ',fontsize=fs)
#plt.xlabel('latitude',fontsize=fs)
ax.set_xlim(right=100)
ax.set_aspect(aspect=1.9)
plt.savefig(savfile+str(variable)+'_alt_VA.png',dpi=300,bbox_inches='tight',format='png')
plt.show()

# pti=[]
# for i in range(alt_max_w,alt_min_w):
#     pti.append(np.mean(dats.var['PT'][i,:,:]))
# 
# # search events
# for i in range(len(track_model3['dates'])): 
#     if track_model3['dates'][i].day==17:
#         print(i)