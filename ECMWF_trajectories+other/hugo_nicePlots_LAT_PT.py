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

mymap=plt.cm.get_cmap('hot', 20)#'RdGy'

# hypothèse :
    #que moyenne TP° chaque jour la même ?
    # moyenne zonale
initials='CAN2017'
nametrackmodel3 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_B2.pkl'
savfile='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/'
variable='PT'

track_model3 = pickle.load(open(nametrackmodel3,'rb'))
date00=track_model3['dates'][0] #choix de l'indice pour l'image de fond
altrangemin=9
altrangemax=25
lnwd=2

if os.path.isfile(savfile+'dats00_PT_'+str(0)):
    with gzip.open(savfile+'dats00_PT_'+str(0),'rb') as f:
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
    PTmean = np.mean(dat.var[variable],axis=2)
    dat.var['PTmean'] = PTmean[:,:,np.newaxis]
    dats00 = dat.extract(varss='All',lonRange=(0,359),latRange=(22,85))
        
    ## Enregistrement des données à réutiliser pour ce plot-ci : 
    with gzip.open(savfile+'dats00_PT_'+str(0),'wb') as f:
        pickle.dump(dats00,f)

## Calcul des isentropes :
lon=track_model3['lons'][0]
ialtrangemin=np.where(dats00.attr['zscale']<=altrangemin)[0][0]
ialtrangemax=np.where(dats00.attr['zscale']<=altrangemax)[0][0]
levs=(ialtrangemax,ialtrangemin)
txt='Meridional view of the zonal mean potential temprature (K°)'
log=False
clim=(None,None)
cmap=mymap
scale=1
try:
    pos=np.where(dats00.attr['lons']>=lon)[0][0]
    print('pos',pos)
except:
    print('lon out of range')
fig2 = plt.figure(figsize=(11,4))
fs = 16
ax2 = fig2.add_subplot(111)
if levs[0]==None: l1=29
else: l1 = levs[0]
if levs[1]==None: l2=115
else: l2 = levs[1]
lats = np.arange(dats00.attr['lats'][0]-0.5*dats00.attr['dla'],dats00.attr['lats'][-1]+dats00.attr['dla'],dats00.attr['dla'])
try:
    zz1 = 0.5*(dats00.var['Z'][l1-1:l2+1,:, pos] + dats00.var['Z'][l1:l2+2,:,pos])/1000
    zz = np.empty((zz1.shape[0],zz1.shape[1]+1))
    zz[:,1:-1] = 0.5*(zz1[:,1:]+zz1[:,:-1])
    zz[:,0] = zz1[:,0]
    zz[:,-1] = zz1[:,-1]
    print(zz.shape,len(lats))
    iax2=ax2.pcolormesh(lats,zz,scale*dats00.var[variable+'mean'][l1:l2+1,:,pos],
                    vmin=clim[0],vmax=clim[1],cmap=cmap)
    plt.ylabel('altitude (km)',fontsize=fs)
    print('USE Z')
except(KeyError):
    iax2=ax2.pcolormesh(lats,dats00.attr['zscale_i'][l1:l2+2],scale*dats00.var[variable+'mean'][l1:l2+1,:, pos],
                vmin=clim[0],vmax=clim[1],cmap=cmap)

    plt.ylabel('baro altitude (km)',fontsize=fs)
ax2.tick_params(labelsize=16)
plt.xlabel('latitude',fontsize=fs)

plt.title(txt,fontsize=fs)#+" longitude="+str(lon)+'°'
cbar = fig2.colorbar(iax2)
cbar.ax.tick_params(labelsize=fs)

## Idée : TP pour niveaux modèles alt -> approcher alt de la trajctoire à ceux des niveaux modèles -> convertir lvl modèle vers TP
# for i in range(len(track_model3['alts'])):
#     # mean sur longitude nécessaire ou juste prendre la valeur à un point donné ? 
#     #lon_min_w=track_model3['lons'][i]-fenetre_lon/2-1
#     #lon_min_wi=np.where(dats00.attr['lons']>=lon_min_w)[0][0]
#     #lon_max_w=track_model3['lons'][i]+fenetre_lon/2+1
#     #lon_max_wi=np.where(dats00.attr['lons']>=lon_max_w)[0][0]
#     lon_wi=np.where(dats00.attr['lons']>=track_model3['lons'][i])[0][0]
#     track_model3['alts'][i]=dats00.var['PT'][np.where(dats00.attr['zscale']<=track_model3['alts'][i])[0][0],np.where(dats00.attr['lats']==np.int(np.round(track_model3['lats'][i])))[0][0],lon_wi]
    
colori=len(track_model3['lats'])
jet3= plt.get_cmap('binary')
colors32 = iter(jet3(np.linspace(0,1,colori)))    

for i in range(0,len(track_model3['lats'])):
    #ax2.plot([track_model3['lats'][i-1],track_model3['lats'][i]],[track_model3['alts'][i-1],track_model3['alts'][i]],linewidth=lnwd,color=next(colors32))
    ax2.scatter(track_model3['lats'][i],track_model3['alts'][i],linewidth=lnwd,color=next(colors32))
plt.savefig(savfile+str(variable)+'_meridional.png',bbox_inches='tight',dpi=300,format='png')
plt.show()  