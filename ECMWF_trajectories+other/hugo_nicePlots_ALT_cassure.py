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
from mpl_toolkits.axes_grid1 import Grid 

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
nametrackmodel3 = '/home/hlestrel/stageM1-vortex/Vortex-track_CAN2017_CASSURE.pkl'
savfile='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/cassure/'
nameexport='/data/hlestrel/fig_bulle_'+initials+'/CONTOUR/export'
with gzip.open(nameexport,'rb') as f:
    export=pickle.load(f)
for i in range(len(export)):
    for j in range(len(export[i]["lons"])):
        if export[i]['lons'][j]>180:
            export[i]['lons'][j]-=360
variable='PV'
# fenetre image fond
lonrangemin=-50
lonrangemax=50
latrangemin=17
latrangemax=85
###### UPDATE TRACK
track_model3 = pickle.load(open(nametrackmodel3,'rb'))

#plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
#plt.rcParams['ytick.left'] = plt.rcParams['ytick.labelleft'] = False

props = dict(boxstyle='round', facecolor='white', alpha=0.8)

# ## double des données pour la coupe zonale :
# for itrack in range(len(track_model3['dates'])):
#     date00=track_model3['dates'][itrack]
#     namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
#     if os.path.isfile(namedate00):
#         with gzip.open(namedate00,'rb') as f:
#             dats=pickle.load(f)
#     else:    
#         dat = ECMWF('FULL-EA',date00,exp='VOZ') 
#         dat._get_var('O3')
#         dat._get_var('T')
#         dat._get_var('U') 
#         dat._get_var('V')
#         dat._get_var('VO')
#         dat._get_var('P')
#         dat._get_var('PT')
#         dat._mkpscale() 
#         dat._mkzscale()
#         dat._mkp() 
#         dat._mkthet() 
#         dat._mkpv()
#         PVmean = np.mean(dat.var[variable],axis=2)
#         dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
#         datp=dat.shift2west(300) 
#         dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
#         with gzip.open(namedate00,'wb') as f:
#             pickle.dump(dats,f)
 
projplate = ccrs.PlateCarree(central_longitude=0)
proj=projplate 
cardinal_level=True
clim=(-10,10)
figsize=(9,22)
fs=15
fig = plt.figure(figsize=figsize,constrained_layout=True)#,sharex='col', sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0})
fontsize=fs
cmap=mymap
scale=10**6
aspect=1
sat_H=35785831
xleg,yleg=proj2legend(-45-lonrangemin,80-latrangemin,lonrangemin,lonrangemax,latrangemin,latrangemax)

##
itrack=0
ax0 = plt.subplot(421,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
if os.path.isfile(namedate00):
    with gzip.open(namedate00,'rb') as f:
        dats=pickle.load(f)
ialts=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
lev=ialts
if len(dats.var[variable+'ano'].shape) == 3:
    if (cardinal_level==False) | (lev > dats.nlev-1):
        clev = np.abs(np.array(dats.attr['levs'])-lev).argmin()
    else:
        clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
iax=ax0.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
xlocs = None
gl = ax0.gridlines(draw_labels=True, xlocs=None,
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
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
gl.rotate_labels=None
if track_model3['dates'][itrack].month==8:
    textstr=str(track_model3['dates'][itrack].day)+' Aug\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
ax0.text(xleg,yleg,textstr, transform=ax0.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax0.label_outer()

##
itrack=1
ax1 = plt.subplot(423,sharex=ax0 ,sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
if os.path.isfile(namedate00):
    with gzip.open(namedate00,'rb') as f:
        dats=pickle.load(f)
ialts=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
lev=ialts
if len(dats.var[variable+'ano'].shape) == 3:
    if (cardinal_level==False) | (lev > dats.nlev-1):
        clev = np.abs(np.array(dats.attr['levs'])-lev).argmin()
    else:
        clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
ax1.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
xlocs = None
gl = ax1.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax1.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax1.coastlines('50m')
gl.top_labels = False
gl.bottom_labels = False
gl.right_labels = False
gl.left_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
gl.rotate_labels=None
if track_model3['dates'][itrack].month==8:
    textstr=str(track_model3['dates'][itrack].day)+' Aug\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
#plt.setp(ax1.get_xticklabels(), visible=False)
ax1.text(xleg,yleg,textstr, transform=ax1.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax1.label_outer()

#
itrack=2
iexp=0
ax2 = plt.subplot(425,sharex=ax0, sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
if os.path.isfile(namedate00):
    with gzip.open(namedate00,'rb') as f:
        dats=pickle.load(f)
ialts=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
lev=ialts
if len(dats.var[variable+'ano'].shape) == 3:
    if (cardinal_level==False) | (lev > dats.nlev-1):
        clev = np.abs(np.array(dats.attr['levs'])-lev).argmin()
    else:
        clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
ax2.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
ax2.plot(export[iexp]['lons'][np.where(export[iexp]['lats']<=latrangemax)[0][0]:np.where(export[iexp]['lats']<=latrangemin)[0][0]],export[iexp]['lats'][np.where(export[iexp]['lats']<=latrangemax)[0][0]:np.where(export[iexp]['lats']<=latrangemin)[0][0]],'k',label='CALIOP orbits')#,linewidth=1,linestyle='-')
xlocs = None
gl = ax2.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax2.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax2.coastlines('50m')
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
gl.rotate_labels=None
if track_model3['dates'][itrack].month==8:
    textstr=str(track_model3['dates'][itrack].day)+' Aug\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
ax2.text(xleg,yleg,textstr, transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax2.label_outer()

#
itrack=3
iexp+=1
ax3 = plt.subplot(427,sharex=ax0, sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
if os.path.isfile(namedate00):
    with gzip.open(namedate00,'rb') as f:
        dats=pickle.load(f)
ialts=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
lev=ialts
if len(dats.var[variable+'ano'].shape) == 3:
    if (cardinal_level==False) | (lev > dats.nlev-1):
        clev = np.abs(np.array(dats.attr['levs'])-lev).argmin()
    else:
        clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
ax3.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
ax3.plot(export[iexp]['lons'][np.where(export[iexp]['lats']<=latrangemax)[0][0]:np.where(export[iexp]['lats']<=latrangemin)[0][0]],export[iexp]['lats'][np.where(export[iexp]['lats']<=latrangemax)[0][0]:np.where(export[iexp]['lats']<=latrangemin)[0][0]],'k')
xlocs = None
gl = ax3.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax3.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax3.coastlines('50m')
gl.top_labels = False
gl.right_labels = False
gl.left_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
gl.rotate_labels=None
if track_model3['dates'][itrack].month==8:
    textstr=str(track_model3['dates'][itrack].day)+' Aug\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
ax3.text(xleg,yleg,textstr, transform=ax3.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax3.label_outer()

##
itrack=4
iexp+=1
ax4 = plt.subplot(422,sharex=ax0, sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
if os.path.isfile(namedate00):
    with gzip.open(namedate00,'rb') as f:
        dats=pickle.load(f)
ialts=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
lev=ialts
if len(dats.var[variable+'ano'].shape) == 3:
    if (cardinal_level==False) | (lev > dats.nlev-1):
        clev = np.abs(np.array(dats.attr['levs'])-lev).argmin()
    else:
        clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
ax4.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
ax4.plot(export[iexp]['lons'][np.where(export[iexp]['lats']<=latrangemax)[0][0]:np.where(export[iexp]['lats']<=latrangemin)[0][0]],export[iexp]['lats'][np.where(export[iexp]['lats']<=latrangemax)[0][0]:np.where(export[iexp]['lats']<=latrangemin)[0][0]],'k')
xlocs = None
gl = ax4.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax4.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax4.coastlines('50m')
gl.top_labels = False
gl.right_labels = False
gl.left_labels = False
gl.bottom_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
if track_model3['dates'][itrack].month==8:
    textstr=str(track_model3['dates'][itrack].day)+' Aug\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
ax4.text(xleg,yleg,textstr, transform=ax4.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax4.label_outer()

##
itrack=5
iexp+=1
ax5 = plt.subplot(424,sharex=ax0, sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
if os.path.isfile(namedate00):
    with gzip.open(namedate00,'rb') as f:
        dats=pickle.load(f)
ialts=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
lev=ialts
if len(dats.var[variable+'ano'].shape) == 3:
    if (cardinal_level==False) | (lev > dats.nlev-1):
        clev = np.abs(np.array(dats.attr['levs'])-lev).argmin()
    else:
        clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
ax5.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
ax5.plot(export[iexp]['lons'][np.where(export[iexp]['lats']<=latrangemax)[0][0]:np.where(export[iexp]['lats']<=latrangemin)[0][0]],export[iexp]['lats'][np.where(export[iexp]['lats']<=latrangemax)[0][0]:np.where(export[iexp]['lats']<=latrangemin)[0][0]],'k')
xlocs = None
gl = ax5.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax5.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax5.coastlines('50m')
gl.top_labels = False
gl.right_labels = True
gl.left_labels = False
gl.bottom_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
if track_model3['dates'][itrack].month==8:
    textstr=str(track_model3['dates'][itrack].day)+' Aug\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
ax5.text(xleg,yleg,textstr, transform=ax5.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax5.label_outer()

##
itrack=6
ax6 = plt.subplot(426,sharex=ax0, sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
if os.path.isfile(namedate00):
    with gzip.open(namedate00,'rb') as f:
        dats=pickle.load(f)
ialts=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
lev=ialts
if len(dats.var[variable+'ano'].shape) == 3:
    if (cardinal_level==False) | (lev > dats.nlev-1):
        clev = np.abs(np.array(dats.attr['levs'])-lev).argmin()
    else:
        clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
ax6.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
xlocs = None
gl = ax6.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax6.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax6.coastlines('50m')
gl.top_labels = False
gl.right_labels = False
gl.left_labels = False
gl.bottom_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
if track_model3['dates'][itrack].month==8:
    textstr=str(track_model3['dates'][itrack].day)+' Aug\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
ax6.text(xleg,yleg,textstr, transform=ax6.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax6.label_outer()

##
itrack=8
iexp+=1
ax7 = plt.subplot(428,sharex=ax0, sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i='+str(itrack)
if os.path.isfile(namedate00):
    with gzip.open(namedate00,'rb') as f:
        dats=pickle.load(f)
ialts=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
lev=ialts
if len(dats.var[variable+'ano'].shape) == 3:
    if (cardinal_level==False) | (lev > dats.nlev-1):
        clev = np.abs(np.array(dats.attr['levs'])-lev).argmin()
    else:
        clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
ax7.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
ax7.plot(export[iexp]['lons'][np.where(export[iexp]['lats']<=latrangemax)[0][0]:np.where(export[iexp]['lats']<=latrangemin)[0][0]],export[iexp]['lats'][np.where(export[iexp]['lats']<=latrangemax)[0][0]:np.where(export[iexp]['lats']<=latrangemin)[0][0]],'k')
xlocs = None
gl = ax7.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax7.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax7.coastlines('50m')
gl.top_labels = False
gl.right_labels = True
gl.left_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
if track_model3['dates'][itrack].month==8:
    textstr=str(track_model3['dates'][itrack].day)+' Aug\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept\nZ = '+str(round(track_model3['alts'][itrack],2))+' km'
ax7.text(xleg,yleg,textstr, transform=ax7.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax7.label_outer()

fig.suptitle('Anomaly of potential vorticity ($10^{-6}$ m$^{2}$ s$^{-1}$ K kg$^{-1}$)',fontsize=fs)
pos_cax = fig.add_axes([1.05, 0.15, 0.05, 0.3])
fig.legend(loc='upper right',fontsize=fs)#, bbox_to_anchor=(0.065, 0.93))
cbar=fig.colorbar(iax,cax=pos_cax)
cbar.ax.tick_params(labelsize=fs) 
plt.subplots_adjust(hspace=0.1,wspace=0.02)
fig.add_subplot(111, frame_on=False)
plt.tick_params(labelcolor="none", bottom=False, left=False)
plt.xlabel('longitude',fontsize=fs)
#plt.tight_layout(h_pad=0,w_pad=0)
fig.savefig(savfile+str(variable)+'_CASSURE_alt.png',dpi=300,bbox_inches='tight',format='png')
plt.show()