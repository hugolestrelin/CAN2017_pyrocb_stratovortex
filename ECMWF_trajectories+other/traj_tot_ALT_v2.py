#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Shows the model from the data of the chosen sequences of CALIOP sections in function of the track dictionary (trackname), from the directory defined in diroutput.
This opration is realised in a separated serveur, the tracks have therefore to be previously copied from a serveur to another. 

input : track dictionary of chosen sequences (diroutput ; trackname)
output : plots of the model in diroutput

@author: Hugo Lestrelin

"""
import sys,os  
import socket      
if 'ciclad' in socket.gethostname():
    import sys,os
    sys.path.append('/home/hlestrel/pylib')                
from datetime import datetime,timedelta
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

initials='CAN2017'
# ERA5 track
savfile='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/fig_traj/'
[trac_PV1,track_model0] = pickle.load(open(savfile+'Vortex0-track.pkl','rb'))
[trac_PV4,track_model1] = pickle.load(open(savfile+'VortexA-track.pkl','rb'))
[trac_PV2,track_model2] = pickle.load(open(savfile+'VortexB1-track.pkl','rb'))
[trac_PV3,track_model3] = pickle.load(open(savfile+'VortexB2-track.pkl','rb'))

itera=1
lnwd=4
s=14

idx0=len(track_model0['dates'])
idx1=len(track_model1['dates'])
idx2=len(track_model2['dates'])
idx3=len(track_model3['dates'])
for i in range(idx0):
    if track_model0['lons'][i]>180:
        track_model0['lons'][i]-=360
for i in range(idx1):
    if track_model1['lons'][i]>180:
        track_model1['lons'][i]-=360
for i in range(idx2):
    if track_model2['lons'][i]>180:
        track_model2['lons'][i]-=360
for i in range(idx3):
    if track_model3['lons'][i]>180:
        track_model3['lons'][i]-=360


colordebut=0.3
colorfin=1
colori0=int(round(len(track_model0['z'])/itera))
colori1=int(round(len(track_model1['z'])/itera))
colori2=int(round(len(track_model2['z'])/itera))
colori3=int(round(len(track_model3['z'])/itera))
jet0= plt.get_cmap('Wistia') # Greys ; winter
colors0 = iter(jet0(np.linspace(colordebut,colorfin,colori0))) 
colors00 = iter(jet0(np.linspace(colordebut,colorfin,colori0))) 
jet1= plt.get_cmap('Reds') #Reds ; spring
colors1 = iter(jet1(np.linspace(colordebut,colorfin,colori1))) 
colors11 = iter(jet1(np.linspace(colordebut,colorfin,colori1))) 
jet2= plt.get_cmap('Blues') #Blues ; summer
colors2 = iter(jet2(np.linspace(colordebut,colorfin,colori2)))
colors22 = iter(jet2(np.linspace(colordebut,colorfin,colori2)))
jet3= plt.get_cmap('Greens') #Greens ; autumn
colors3 = iter(jet3(np.linspace(colordebut,colorfin,colori3)))
colors33 = iter(jet3(np.linspace(colordebut,colorfin,colori3)))

## COUPE LONGITUDINALE
fig2 = plt.figure(figsize=(11,4))
fs = 18
ax2 = fig2.add_subplot(111)

plt.ylabel('Potential Temperature (K)',fontsize=fs)
ax2.tick_params(labelsize=fs)
plt.xlabel('longitude',fontsize=fs)

idx=idx0
dates=track_model0['dates']
lons1=track_model0['lons']
alts1=track_model0['pt']
for i in range(itera,idx,itera):
    c=next(colors0)
    if  (dates[i].month==9 and 3<dates[i].day<16):
        ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
    else:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax2.plot([180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
        elif (dates[i].month==9 and dates[i].day==1 and dates[i].hour==0):
            ax2.plot([lons1[len(lons1)-1-itera],lons1[len(lons1)-1]],[alts1[len(lons1)-1-itera],alts1[len(lons1)-1]],linewidth=lnwd,color=c,label='Vortex 0')
        elif (lons1[i-itera]>150 and lons1[i]<-150):
            ax2.plot([-180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
        else:
            ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)

idx=idx1
dates=track_model1['dates']
lons1=track_model1['lons']
alts1=track_model1['pt']
for i in range(itera,idx,itera):
    c=next(colors1)
    if  (dates[i].month==9 and 3<dates[i].day<16):
        ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
    else:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax2.plot([180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
        elif (dates[i].month==9 and dates[i].day==16 and dates[i].hour==0):
            ax2.plot([lons1[len(lons1)-1-itera],lons1[len(lons1)-1]],[alts1[len(lons1)-1-itera],alts1[len(lons1)-1]],linewidth=lnwd,color=c,label='Vortex A')
        elif (lons1[i-itera]>150 and lons1[i]<-150):
            ax2.plot([-180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
        else:
            ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
            
idx=idx2
dates=track_model2['dates']
lons1=track_model2['lons']
alts1=track_model2['pt']
for i in range(itera,idx,itera):
    c=next(colors2)
    if  (dates[i].month==9 and 3<dates[i].day<16):
        ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
    else:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax2.plot([180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
        elif (dates[i].month==9 and dates[i].day==16 and dates[i].hour==0):
            ax2.plot([lons1[len(lons1)-1-itera],lons1[len(lons1)-1]],[alts1[len(lons1)-1-itera],alts1[len(lons1)-1]],linewidth=lnwd,color=c,label='Vortex B1')
        elif (lons1[i-itera]>150 and lons1[i]<-150):
            ax2.plot([-180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
        else:
            ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)

idx=idx3
dates=track_model3['dates']
lons1=track_model3['lons']
alts1=track_model3['pt']
for i in range(itera,idx,itera):
    c=next(colors3)
    if  (dates[i].month==9 and 3<dates[i].day<16):
        if (dates[i].day==15 and dates[i].hour==0):
            ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c,marker='x',label='CALIOP down',markersize=s)
        else:
            ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
    else:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax2.plot([180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
        elif (dates[i].month==9 and dates[i].day==16 and dates[i].hour==0):
            ax2.plot([lons1[len(lons1)-1-itera],lons1[len(lons1)-1]],[alts1[len(lons1)-1-itera],alts1[len(lons1)-1]],linewidth=lnwd,color=c,label='Vortex B2')
        elif (lons1[i-itera]>150 and lons1[i]<-150):
            ax2.plot([-180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
        else:
            ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
            
ax2.text(-0.15,1.05,'a)',transform=ax2.transAxes,fontsize=16)
plt.xticks([-180,-120,-60,0,60,120,180],['180°W','120°W','60°W','0°','60°E','120°E','180°E'])
ax2.set_xlim(-180,180)
plt.legend(bbox_to_anchor=(0.745, 0.527), loc='upper left', borderaxespad=0.,fontsize=15)
plt.savefig(savfile+'_lon_PT_'+initials+'.png',bbox_inches='tight',dpi=300,format='png')

## PLOT HORIZONTAL
scale=10**6
projplate = ccrs.PlateCarree(central_longitude=0)
proj=projplate
fig1 = plt.figure(figsize=(11,4))
fig1.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
ax = plt.axes(projection = proj)
gl = ax.gridlines(draw_labels=True,
                linewidth=2, color='gray', alpha=0.5, linestyle='--')
ax.add_feature(feature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none'))
ax.coastlines('50m')
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': fs}
gl.ylabel_style = {'size': fs}
ax.set_aspect(aspect=1.7)

idx=idx0
dates=track_model0['dates']
lons1=track_model0['lons']
lats1=track_model0['lats']
for i in range(itera,idx,itera):
    c=next(colors00)
    if  dates[i].month==9 and 3<dates[i].day<16:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax.plot([180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
        elif (lons1[i-itera]>150 and lons1[i])<-150:
            ax.plot([-180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
        else:
            ax.plot([lons1[i-itera],lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
    else:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax.plot([180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
        elif (lons1[i-itera]>150 and lons1[i])<-150:
            ax.plot([-180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
        else:
            ax.plot([lons1[i-itera],lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)

idx=idx1
dates=track_model1['dates']
lons1=track_model1['lons']
lats1=track_model1['lats']
for i in range(itera,idx,itera):
    c=next(colors11)
    if  dates[i].month==9 and 3<dates[i].day<16:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax.plot([180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
        elif (lons1[i-itera]>150 and lons1[i])<-150:
            ax.plot([-180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
        else:
            ax.plot([lons1[i-itera],lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
    else:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax.plot([180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
        elif (lons1[i-itera]>150 and lons1[i])<-150:
            ax.plot([-180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
        else:
            ax.plot([lons1[i-itera],lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
            
            
idx=idx2
dates=track_model2['dates']
lons1=track_model2['lons']
lats1=track_model2['lats']
for i in range(itera,idx,itera):
    c=next(colors22)
    if  dates[i].month==9 and 3<dates[i].day<16:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax.plot([180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
        elif (lons1[i-itera]>150 and lons1[i])<-150:
            ax.plot([-180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
        else:
            ax.plot([lons1[i-itera],lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
    else:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax.plot([180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
        elif (lons1[i-itera]>150 and lons1[i])<-150:
            ax.plot([-180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
        else:
            ax.plot([lons1[i-itera],lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)

idx=idx3
dates=track_model3['dates']
lons1=track_model3['lons']
lats1=track_model3['lats']
for i in range(itera,idx,itera):
    c=next(colors33)
    if  dates[i].month==9 and 3<dates[i].day<16:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax.plot([180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
        elif (lons1[i-itera]>150 and lons1[i])<-150:
            ax.plot([-180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
        else:
            ax.plot([lons1[i-itera],lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x',markersize=s)
    else:
        if (lons1[i-itera]-lons1[i]<-100): 
            ax.plot([180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
        elif (lons1[i-itera]>150 and lons1[i])<-150:
            ax.plot([-180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
        else:
            ax.plot([lons1[i-itera],lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)

ax.text(-0.15,1.05,'a)',transform=ax.transAxes,fontsize=16)

plt.savefig(savfile+'_horiz_PT_'+initials+'.png',dpi=300,bbox_inches='tight',format='png')
#plt.show()
