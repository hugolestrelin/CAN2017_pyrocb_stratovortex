#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

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
nametrackmodel3 = '/home/hlestrel/stageM1-vortex/CAN2017/Vortex-track_CAN2017_CASSURE.pkl'
savfile='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/cassure/'
nameexport='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/cassure/export-cassure'
with gzip.open(nameexport,'rb') as f:
    export=pickle.load(f)
with gzip.open('/data/hlestrel/fig_bulle_'+initials+'/fig_publi/cassure/Ntrack-cassure','rb') as f:
    CALIOPtrack=pickle.load(f)
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
latrangemaxlats=69
lvlcon=np.arange(-5,-2,3)
###### UPDATE TRACK
track_model3 = pickle.load(open(nametrackmodel3,'rb'))

props = dict(boxstyle='round', facecolor='white', alpha=0.8)

lnwd=3
projplate = ccrs.PlateCarree(central_longitude=0)
proj=projplate 
cardinal_level=True
clim=(-5,0)
figsize=(9,22)
fs=15
fig = plt.figure(figsize=figsize,constrained_layout=True)#,sharex='col', sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0})
fontsize=fs
cmap1=plt.cm.get_cmap('Blues', 10)#mymap
cmap=cmap1.reversed()
jet2=plt.cm.get_cmap('Wistia',10)#mymap
cmap2=jet2.reversed()
#cmap2=cmap2.reversed()
pair=False
if pair:
    cmaplist = [cmap(i) for i in range(cmap.N)]
    for i in range(0,len(cmaplist),2):
        a=cmaplist[i]
        cmaplist[i]=cmaplist[i+1]
        cmaplist[i+1]=a
    cmap = colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
scale=10**6
aspect=1
sat_H=35785831
xleg,yleg=proj2legend(-45-lonrangemin,80-latrangemin,lonrangemin,lonrangemax,latrangemin,latrangemax)

## 31/08
itrack=0
ax0 = plt.subplot(528,projection = proj)
date00=datetime(2017, 8, 31, 18, 0)
track_model3['dates'][itrack]=date00#datetime(2017, 8, 31, 3, 0)
track_model3['alts'][itrack]=18.5
track_model3['lats'][itrack]=54.0
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_PT_i='+str(itrack)
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
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    datp=dat.shift2west(300) 
    dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
    with gzip.open(namedate00,'wb') as f:
        pickle.dump(dats,f)
dats.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
pos=np.where(dats.attr['lats']>=track_model3['lats'][itrack])[0][0]
lev=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
PTref=465#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
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
##bis
PTref2=460#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref2)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
#ax0.imshow(scale*buf, transform=projplate,cmap=cmap2,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect,alpha=0.5)
#ax0.contour(scale*buf[:,],levels=lvlcon, transform=projplate,cmap=cmap2,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],linewidths=3,alpha=1)
##bis
for j in range(len(export)):
    if CALIOPtrack[j]['date'].day==date00.day:
        if j==14:
            continue
        elif j==15:
            if export[j]['lats'][0]>export[j]['lats'][1]:
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']<=67)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=67)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            elif export[j]['lats'][0]<export[j]['lats'][1]:
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=67)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=67)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
        elif j==37:
            if export[j]['lats'][0]>export[j]['lats'][1]:
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'darkviolet',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            elif export[j]['lats'][0]<export[j]['lats'][1]:
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'darkviolet',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
        else:
            if export[j]['lats'][0]>export[j]['lats'][1]:
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            elif export[j]['lats'][0]<export[j]['lats'][1]:
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax0.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
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
gl.left_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
gl.rotate_labels=None
if track_model3['dates'][itrack].month==8:
    textstr='i)\n'+str(track_model3['dates'][itrack].day)+' Aug '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
elif track_model3['dates'][itrack].month==9:
    textstr='i)\n'+str(track_model3['dates'][itrack].day)+' Sept '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
ax0.text(xleg,yleg,textstr, transform=ax0.transAxes, fontsize=12, verticalalignment='top', bbox=props)
#ax0.text(-0.15,1.05,'h)',transform=ax0.transAxes,fontsize=16)
ax0.label_outer()

## 22/08
itrack=1
ax1 = plt.subplot(521,projection = proj)#,sharex=ax0 ,sharey=ax0)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i_PT='+str(itrack)
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
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    datp=dat.shift2west(300) 
    dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
    with gzip.open(namedate00,'wb') as f:
        pickle.dump(dats,f)
dats.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
pos=np.where(dats.attr['lats']>=track_model3['lats'][itrack])[0][0]
lev=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
PTref=420#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
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
##bis
PTref2=425#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref2)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
#ax1.imshow(scale*buf, transform=projplate,cmap=cmap2,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect,alpha=0.5)
#ax1.contour(scale*buf,levels=lvlcon, transform=projplate,cmap=cmap2,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],linewidths=3,alpha=1,label='VA - '+str(PTref2))
##bis
for j in range(len(export)):
    if CALIOPtrack[j]['date'].day==date00.day:
        if export[j]['lats'][0]>export[j]['lats'][1]:
            ax1.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')
            ax1.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
        elif export[j]['lats'][0]<export[j]['lats'][1]:
            ax1.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
            ax1.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
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
#gl.left_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
gl.rotate_labels=None
if track_model3['dates'][itrack].month==8:
    textstr='a)\n'+str(track_model3['dates'][itrack].day)+' Aug '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
#plt.setp(ax1.get_xticklabels(), visible=False)
ax1.text(xleg,yleg,textstr, transform=ax1.transAxes, fontsize=12, verticalalignment='top', bbox=props)
#ax1.text(-0.15,1.05,'a)',transform=ax1.transAxes,fontsize=16)
ax1.label_outer()

## 23/08
itrack=2
ax2 = plt.subplot(523,projection = proj)#,sharex=ax0 ,sharey=ax0)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i_PT='+str(itrack)
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
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    datp=dat.shift2west(300) 
    dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
    with gzip.open(namedate00,'wb') as f:
        pickle.dump(dats,f)
dats.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
pos=np.where(dats.attr['lats']>=track_model3['lats'][itrack])[0][0]
lev=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
PTref=420#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
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
##bis
PTref2=435#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref2)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
#ax2.imshow(scale*buf, transform=projplate,cmap=cmap2,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect,alpha=0.5)
ax2.contour(scale*buf[:,:50],levels=lvlcon, transform=projplate,cmap=cmap2,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][50]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],linewidths=3,alpha=1,label='VA - '+str(PTref2))
##bis
## ORBITES MULTIPLES
for j in range(len(export)):
    if CALIOPtrack[j]['date'].day==date00.day:
        if export[j]['lats'][0]>export[j]['lats'][1]:
            ax2.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')
            ax2.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)
        elif export[j]['lats'][0]<export[j]['lats'][1]:
            ax2.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
            ax2.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
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
    textstr='b)\n'+str(track_model3['dates'][itrack].day)+' Aug '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K\n('+str(PTref2)+')'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K\n('+str(PTref2)+')'
ax2.text(xleg,yleg,textstr, transform=ax2.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax2.label_outer()

## 25/08
itrack=3
ax3 = plt.subplot(525,sharex=ax0, sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i_PT='+str(itrack)
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
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    datp=dat.shift2west(300) 
    dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
    with gzip.open(namedate00,'wb') as f:
        pickle.dump(dats,f)
dats.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
pos=np.where(dats.attr['lats']>=track_model3['lats'][itrack])[0][0]
lev=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
PTref=425#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
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
##bis
PTref2=455#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref2)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
#ax3.imshow(scale*buf, transform=projplate,cmap=cmap2,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect,alpha=0.5)
ax3.contour(scale*buf[:,:60],levels=lvlcon, transform=projplate,cmap=cmap2,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][60]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],linewidths=3,alpha=1,label='VA - '+str(PTref2))
##bis
for j in range(len(export)):
    if CALIOPtrack[j]['date'].day==date00.day:
        if j==26:
            if export[j]['lats'][0]>export[j]['lats'][1]:
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'darviolet',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            elif export[j]['lats'][0]<export[j]['lats'][1]:
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'darkviolet',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
        elif j==23:
            if export[j]['lats'][0]>export[j]['lats'][1]:
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'limegreen',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            elif export[j]['lats'][0]<export[j]['lats'][1]:
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'limegreen',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
        else:
            if export[j]['lats'][0]>export[j]['lats'][1]:
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            elif export[j]['lats'][0]<export[j]['lats'][1]:
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax3.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
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
#gl.left_labels = False
gl.bottom_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
gl.rotate_labels=None
if track_model3['dates'][itrack].month==8:
    textstr='c)\n'+str(track_model3['dates'][itrack].day)+' Aug '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K\n('+str(PTref2)+')'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K\n('+str(PTref2)+')'
ax3.text(xleg,yleg,textstr, transform=ax3.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax3.label_outer()

## 26/08
itrack=4
ax4 = plt.subplot(527,sharex=ax0, sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i_PT='+str(itrack)
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
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    datp=dat.shift2west(300) 
    dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
    with gzip.open(namedate00,'wb') as f:
        pickle.dump(dats,f)
dats.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
pos=np.where(dats.attr['lats']>=track_model3['lats'][itrack])[0][0]
lev=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
PTref=435#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
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
##bis
PTref2=460#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref2)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
#ax4.imshow(scale*buf, transform=projplate,cmap=cmap2,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect,alpha=0.5)
ax4.contour(scale*buf,levels=lvlcon, transform=projplate,cmap=cmap2,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],linewidths=3,alpha=1,label='VA - '+str(PTref2))
##bis
for j in range(len(export)):
    if CALIOPtrack[j]['date'].day==date00.day:
        if export[j]['lats'][0]>export[j]['lats'][1]:
            ax4.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')
            ax4.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
        elif export[j]['lats'][0]<export[j]['lats'][1]:
            ax4.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
            ax4.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
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
#gl.left_labels = False
gl.bottom_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
if track_model3['dates'][itrack].month==8:
    textstr='d)\n'+str(track_model3['dates'][itrack].day)+' Aug '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K\n('+str(PTref2)+')'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K\n('+str(PTref2)+')'
ax4.text(xleg,yleg,textstr, transform=ax4.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax4.label_outer()

## 27/08
itrack=9
ax9 = plt.subplot(529,sharex=ax0, sharey=ax0,projection = proj)
track_model3['dates'].append(datetime(2017, 8, 27, 9, 0))
track_model3['alts'].append(17.3)
track_model3['lats'].append(47.0)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i_PT='+str(itrack)
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
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    datp=dat.shift2west(300) 
    dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
    with gzip.open(namedate00,'wb') as f:
        pickle.dump(dats,f)
dats.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
pos=np.where(dats.attr['lats']>=track_model3['lats'][itrack])[0][0]
lev=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
PTref=440#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
ax9.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
##bis
PTref2=460#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref2)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
#ax9.imshow(scale*buf, transform=projplate,cmap=cmap2,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect,alpha=0.5)
ax9.contour(scale*buf,levels=lvlcon, transform=projplate,cmap=cmap2,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],linewidths=3,alpha=1,label='VA - '+str(PTref2))
##bis
for j in range(len(export)):
    if CALIOPtrack[j]['date'].day==date00.day:
        if j==6:
            if export[j]['lats'][0]>export[j]['lats'][1]:
                ax9.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'darkviolet',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax9.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            elif export[j]['lats'][0]<export[j]['lats'][1]:
                ax9.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'darkviolet',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax9.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
        else:
            if export[j]['lats'][0]>export[j]['lats'][1]:
                ax9.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax9.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            elif export[j]['lats'][0]<export[j]['lats'][1]:
                ax9.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax9.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
gl = ax9.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax9.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax9.coastlines('50m')
gl.top_labels = False
gl.right_labels = False
#gl.left_labels = False
#gl.bottom_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
if track_model3['dates'][itrack].month==8:
    textstr='e)\n'+str(track_model3['dates'][itrack].day)+' Aug '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K\n('+str(PTref2)+')'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K\n('+str(PTref2)+')'
ax9.text(xleg,yleg,textstr, transform=ax9.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax9.label_outer()

## 28/08
itrack=5
ax5 = plt.subplot(522,sharex=ax0, sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i_PT='+str(itrack)
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
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    datp=dat.shift2west(300) 
    dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
    with gzip.open(namedate00,'wb') as f:
        pickle.dump(dats,f)
dats.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
pos=np.where(dats.attr['lats']>=track_model3['lats'][itrack])[0][0]
lev=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
PTref=450#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
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
##bis
PTref2=465#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref2)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
#ax5.imshow(scale*buf, transform=projplate,cmap=cmap2,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect,alpha=0.5)
#ax5.contour(scale*buf,levels=lvlcon, transform=projplate,cmap=cmap2,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],linewidths=3,alpha=1,label='VA - '+str(PTref2))
##bis
for j in range(len(export)):
    if CALIOPtrack[j]['date'].day==date00.day:
        if export[j]['lats'][0]>export[j]['lats'][1]:
            ax5.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')
            ax5.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
        elif export[j]['lats'][0]<export[j]['lats'][1]:
            ax5.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
            ax5.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
gl = ax5.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax5.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax5.coastlines('50m')
gl.top_labels = False
gl.right_labels = False
gl.left_labels = False
gl.bottom_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
if track_model3['dates'][itrack].month==8:
    textstr='f)\n'+str(track_model3['dates'][itrack].day)+' Aug '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
ax5.text(xleg,yleg,textstr, transform=ax5.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax5.label_outer()

## 29/08
itrack=10
ax10 = plt.subplot(524,sharex=ax0, sharey=ax0,projection = proj)
track_model3['dates'].append(datetime(2017, 8, 29, 3, 0))
track_model3['alts'].append(18.0)
track_model3['lats'].append(46.0)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i_PT='+str(itrack)
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
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    datp=dat.shift2west(300) 
    dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
    with gzip.open(namedate00,'wb') as f:
        pickle.dump(dats,f)
dats.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
pos=np.where(dats.attr['lats']>=track_model3['lats'][itrack])[0][0]
lev=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
PTref=455#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
ax10.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
##bis
PTref2=465#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref2)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
#ax10.imshow(scale*buf, transform=projplate,cmap=cmap2,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect,alpha=0.5)
#ax10.contour(scale*buf,levels=lvlcon, transform=projplate,cmap=cmap2,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],linewidths=3,alpha=1,label='VA - '+str(PTref2))
##bis
for j in range(len(export)):
    if CALIOPtrack[j]['date'].day==date00.day:
        if j==10 or j==33:
            if export[j]['lats'][0]>export[j]['lats'][1]:
                ax10.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'limegreen',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax10.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            elif export[j]['lats'][0]<export[j]['lats'][1]:
                ax10.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'limegreen',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax10.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
        else:
            if export[j]['lats'][0]>export[j]['lats'][1]:
                ax10.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax10.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            elif export[j]['lats'][0]<export[j]['lats'][1]:
                ax10.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax10.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
gl = ax10.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax10.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax10.coastlines('50m')
gl.top_labels = False
gl.right_labels = False
gl.left_labels = False
gl.bottom_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
if track_model3['dates'][itrack].month==8:
    textstr='g)\n'+str(track_model3['dates'][itrack].day)+' Aug '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
ax10.text(xleg,yleg,textstr, transform=ax10.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax10.label_outer()

## 30/08
itrack=6
ax6 = plt.subplot(526,sharex=ax0, sharey=ax0,projection = proj)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i_PT='+str(itrack)
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
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    datp=dat.shift2west(300) 
    dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
    with gzip.open(namedate00,'wb') as f:
        pickle.dump(dats,f)
dats.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
pos=np.where(dats.attr['lats']>=track_model3['lats'][itrack])[0][0]
lev=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
PTref=470#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
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
##bis
PTref2=460#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref2)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
#ax6.imshow(scale*buf, transform=projplate,cmap=cmap2,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect,alpha=0.5)
#ax6.contour(scale*buf,levels=lvlcon, transform=projplate,cmap=cmap2,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],linewidths=3,alpha=1,label='VA - '+str(PTref2))
##bis
for j in range(len(export)):
    if CALIOPtrack[j]['date'].day==date00.day:
        if export[j]['lats'][0]>export[j]['lats'][1]:
            ax6.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')
            ax6.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
        elif export[j]['lats'][0]<export[j]['lats'][1]:
            if j==19:
                ax6.plot(export[j]['lons'][np.where(export[j]['lats']>=56.5)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=56.5)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
            else:
                ax6.plot(export[j]['lons'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],export[j]['lats'][np.where(export[j]['lats']>=latrangemin)[0][0]:np.where(export[j]['lats']>=latrangemaxlats)[0][0]],'k',linewidth=lnwd)#,linewidth=1,linestyle='-')xlocs = None
                ax6.plot(export[j]['lons'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],export[j]['lats'][np.where(export[j]['lats']>=CALIOPtrack[j]['south'])[0][0]:np.where(export[j]['lats']>=CALIOPtrack[j]['north'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
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
    textstr='h)\n'+str(track_model3['dates'][itrack].day)+' Aug '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
elif track_model3['dates'][itrack].month==9:
    textstr=str(track_model3['dates'][itrack].day)+' Sept '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
ax6.text(xleg,yleg,textstr, transform=ax6.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax6.label_outer()

## 01/09
itrack=8
ax7 = fig.add_subplot(5,2,10,sharex=ax0, sharey=ax0,projection = proj)#plt.subplot(5210,sharex=ax0, sharey=ax0,projection = proj)
track_model3['dates'][itrack]=datetime(2017, 9, 1, 0, 0)
date00=track_model3['dates'][itrack]
namedate00=savfile+'dats_'+str(date00.day)+'_'+str(date00.month)+'_i_PT='+str(itrack)
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
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    dat.var[variable+'ano'] = dat.var[variable] - PVmean[:,:,np.newaxis]
    datp=dat.shift2west(300) 
    dats = datp.extract(varss='All',lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
    with gzip.open(namedate00,'wb') as f:
        pickle.dump(dats,f)
dats.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
pos=np.where(dats.attr['lats']>=track_model3['lats'][itrack])[0][0]
lev=np.where(dats.attr['zscale']<=track_model3['alts'][itrack])[0][0]
PTref=465#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
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
##bis
PTref2=460#dats.var['PTmean'][lev,pos]
for i in range(len(dats.attr['lats'])):
    PTlev=np.where(dats.var['PTmean'][:,i]<=PTref2)[0][0]
    dats.var[variable+'ano'][lev,i,:]=dats.var[variable+'ano'][PTlev,i,:]
if len(dats.var[variable+'ano'].shape) == 3:
    clev = lev
    buf = dats.var[variable+'ano'][clev,:,:]
else:
    buf = dats.var[variable+'ano']
cm_lon =0
ctrl_lon=(dats.attr['lons'][0]+dats.attr['lons'][-1])/2
ctrl_lat=(dats.attr['lats'][0]+dats.attr['lats'][-1])/2
if dats.attr['lons'][-1] > 180: cm_lon=180
projplate = ccrs.PlateCarree(central_longitude=cm_lon)
#ax7.imshow(scale*buf, transform=projplate,cmap=cmap2,clim=clim,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect,alpha=0.5)
#ax7.contour(scale*buf,levels=lvlcon, transform=projplate,cmap=cmap2,extent=[dats.attr['lons'][0]-cm_lon, dats.attr['lons'][-1]-cm_lon,dats.attr['lats'][0], dats.attr['lats'][-1]],linewidths=3,alpha=1,label='VA - '+str(PTref2))
##bis
for j in range(len(export)):
    if CALIOPtrack[j]['date'].day==date00.day:
        if export[j]['lats'][0]>export[j]['lats'][1]:
            if j==18:
                #ax7.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'darkorange',linewidth=lnwd,label='PV anomaly contours at -2 PVU')#,linewidth=1,linestyle='-')
                ax7.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
            else:
                ax7.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'k',label='CALIOP orbits',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax7.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',label='selected occurences',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax7.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'darkviolet',label='detailed orbits \nin Figure 1',linewidth=lnwd)#,linewidth=1,linestyle='-')
                ax7.plot(export[j]['lons'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],export[j]['lats'][np.where(export[j]['lats']<=latrangemaxlats)[0][0]:np.where(export[j]['lats']<=latrangemin)[0][0]],'limegreen',linewidth=lnwd,label='detailed orbits \nin Figure 4')#,linewidth=1,linestyle='-')
                ax7.plot(export[j]['lons'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],export[j]['lats'][np.where(export[j]['lats']<=CALIOPtrack[j]['north'])[0][0]:np.where(export[j]['lats']<=CALIOPtrack[j]['south'])[0][0]],'red',linewidth=lnwd)#,linewidth=1,linestyle='-')
gl = ax7.gridlines(draw_labels=True, xlocs=None,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax7.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax7.coastlines('50m')
gl.top_labels = False
gl.right_labels = False
gl.left_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
if track_model3['dates'][itrack].month==8:
    textstr='j)\n'+str(track_model3['dates'][itrack].day)+' Aug '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
elif track_model3['dates'][itrack].month==9:
    textstr='j)\n'+str(track_model3['dates'][itrack].day)+' Sept '+str(track_model3['dates'][itrack].hour)+' UTC\nPT = '+str(int(PTref))+' K'#\n('+str(PTref2)+')'
ax7.text(xleg,yleg,textstr, transform=ax7.transAxes, fontsize=12, verticalalignment='top', bbox=props)
ax7.label_outer()

fig.suptitle('Anomaly of potential vorticity ($10^{-6}$ m$^{2}$ s$^{-1}$ K kg$^{-1}$)',fontsize=fs)
pos_cax = fig.add_axes([1.05, 0.15, 0.05, 0.3])
fig.legend(loc='upper right',fontsize=fs)#, bbox_to_anchor=(0.065, 0.93))
cbar=fig.colorbar(iax,cax=pos_cax)#orientation="horizontal", pad=0.2)#,cax=pos_cax)
cbar.ax.tick_params(labelsize=fs) 
fig.subplots_adjust(hspace=0.0,wspace=0.02)
fig.add_subplot(111, frame_on=False)
plt.tick_params(labelcolor="none", bottom=False, left=False)
plt.xlabel('longitude',fontsize=fs)
#plt.tight_layout(h_pad=0,w_pad=0)
fig.savefig(savfile+str(variable)+'_CASSURE_PT.png',dpi=300,bbox_inches='tight',format='png')
plt.show()