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
from ECMWF_N import ECMWF 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm
import pickle,gzip
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.dates as mdates
from matplotlib.lines import Line2D

marker_style = dict(linestyle=':', color='0.8', markersize=10,
                    mfc="C0", mec="C0")

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

initials='CAN2017'
savfile='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/fig_traj/'
[trac_PV1,track_model0] = pickle.load(open(savfile+'Vortex0-track.pkl','rb'))
[trac_PV4,track_model1] = pickle.load(open(savfile+'VortexA-track.pkl','rb'))
[trac_PV2,track_model2] = pickle.load(open(savfile+'VortexB1-track.pkl','rb'))
[trac_PV3,track_model3] = pickle.load(open(savfile+'VortexB2-track.pkl','rb'))

itera=6
colordebut=0.3
colorfin=1
colori0=int(round(len(track_model0['z'])/itera))
colori1=int(round(len(track_model1['z'])/itera))
colori2=int(round(len(track_model2['z'])/itera))
colori3=int(round(len(track_model3['z'])/itera))
jet0= plt.get_cmap('Wistia') # Greys ; winter
colors02 = iter(jet0(np.linspace(colordebut,colorfin,colori0))) 
colors0 = iter(jet0(np.linspace(colordebut,colorfin,colori0))) 
colors02d = iter(jet0(np.linspace(0.1,0.15,10))) 
colors0d = iter(jet0(np.linspace(0.1,0.15,10))) 
jet1= plt.get_cmap('Reds') #Reds ; spring
colors12 = iter(jet1(np.linspace(colordebut,colorfin,colori1))) 
colors1 = iter(jet1(np.linspace(colordebut,colorfin,colori1))) 
colors12d = iter(jet1(np.linspace(0.1,0.15,10))) 
colors1d = iter(jet1(np.linspace(0.1,0.15,10))) 
jet2= plt.get_cmap('Blues') #Blues ; summer
colors22 = iter(jet2(np.linspace(colordebut,colorfin,colori2))) 
colors2 = iter(jet2(np.linspace(colordebut,colorfin,colori2)))
colors22d = iter(jet2(np.linspace(0.1,0.15,colori2))) 
colors2d = iter(jet2(np.linspace(0.1,0.15,colori2)))
jet3= plt.get_cmap('Greens') #Greens ; autumn
colors32 = iter(jet3(np.linspace(colordebut,colorfin,colori3))) 
colors3 = iter(jet3(np.linspace(colordebut,colorfin,colori3)))
colors32d = iter(jet3(np.linspace(0.1,0.15,10))) 
colors3d = iter(jet3(np.linspace(0.1,0.15,10)))

fig, ax = plt.subplots(figsize=(11, 4))
plt.ylim(11,40)

ax2=ax.twinx()

lnwd=4
scission=30
ls='-'
mk1='o'
mk2='D'

## Pour la legende :
ax.plot(track_model0['dates'][0],track_model0['z'][0],linewidth=lnwd,color=jet0(np.linspace(colordebut,1,colori0))[int(round(colori0/2))],label='Vortex 0')#color=next(colors12)
ax.plot(track_model1['dates'][0],track_model1['z'][0],linewidth=lnwd,color=jet1(np.linspace(colordebut,1,colori1))[int(round(colori1/2))],label='Vortex A')#color=next(colors12)
ax.plot(track_model2['dates'][0],track_model2['z'][0],linewidth=lnwd,color=jet2(np.linspace(colordebut,1,colori2))[int(round((colori2)/2))],label='Vortex B1')#color=next(colors22)
ax.plot(track_model3['dates'][0],track_model3['z'][0],linewidth=lnwd,color=jet3(np.linspace(colordebut,1,colori3))[int(round((colori3)/2))],label='Vortex B2')#color=next(colors32)  
  
for i in range(itera,len(track_model0['z']),itera):
    c0=next(colors02)
    ax.plot([track_model0['dates'][i-itera],track_model0['dates'][i]],[track_model0['z'][i-itera],track_model0['z'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=c0)
    ax2.plot([track_model0['dates'][i-itera],track_model0['dates'][i]],[track_model0['pt'][i-itera],track_model0['pt'][i]],linewidth=lnwd,marker=mk2,color=next(colors0))#, linestyle='--')
    
ax.plot([track_model0['dates'][i-itera],track_model0['dates'][i]],[track_model0['z'][i-itera],track_model0['z'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=c0,label='Altitude')
ax2.plot([track_model0['dates'][i-itera],track_model0['dates'][i]],[track_model0['pt'][i-itera],track_model0['pt'][i]],linewidth=lnwd,marker=mk2,color=c0,label='Potential temperature')#, linestyle='--')
    
for i in range(itera,len(track_model1['z']),itera):
    ax.plot([track_model1['dates'][i-itera],track_model1['dates'][i]],[track_model1['z'][i-itera],track_model1['z'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=next(colors12))
    ax2.plot([track_model1['dates'][i-itera],track_model1['dates'][i]],[track_model1['pt'][i-itera],track_model1['pt'][i]],linewidth=lnwd,marker=mk2,color=next(colors1))#, linestyle='--')
    
for i in range(itera,len(track_model2['z']),itera):
    ax.plot([track_model2['dates'][i-itera],track_model2['dates'][i]],[track_model2['z'][i-itera],track_model2['z'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=next(colors22))
    ax2.plot([track_model2['dates'][i-itera],track_model2['dates'][i]],[track_model2['pt'][i-itera],track_model2['pt'][i]],linewidth=lnwd,marker=mk2,color=next(colors2))#, linestyle='--')
    
    
for i in range(itera,len(track_model3['z']),itera):
    c=next(colors32)
    ax.plot([track_model3['dates'][i-itera],track_model3['dates'][i]],[track_model3['z'][i-itera],track_model3['z'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=c)
    ax2.plot([track_model3['dates'][i-itera],track_model3['dates'][i]],[track_model3['pt'][i-itera],track_model3['pt'][i]],linewidth=lnwd,marker=mk2,color=next(colors3))#, linestyle='--')
    #ax2.scatter(track_model3['lats'][i],track_model3['alts'][i],linewidth=lnwd,color=next(colors32))
    
#ax.plot(track_model3['dates'][-1],track_model3['z'][-1],linewidth=lnwd,color=jet3(np.linspace(colordebut,1,colori3))[-1],marker='o',label='Altitude')#color=next(colors12)
#ax2.plot(track_model3['dates'][-1],track_model3['pt'][-1],linewidth=lnwd,color=jet3(np.linspace(colordebut,1,colori3))[-1],marker='D',label='Potential temperature')#color=next(colors22)

fs=18
#ax.set(xlabel="Date",ylabel="Altitude (km)") #baro
#ax.xlabel_style = {'fontsize': fs}
#ax.ylabel_style = {'fontsize': fs}
ax.set_xlabel("Date",fontsize=fs) #baro
ax.set_ylabel("Altitude (km)",fontsize=fs) #baro
ax.set_yticks([15,20,25])#,['15','20','25'])

ax2.set_ylabel("Potential Temperature (K)",fontsize=fs) #baro
#ax2.ylabel_style = {'size': fs}
#ax.legend(loc='upper left')
#fig.legend(loc='upper left', bbox_to_anchor=(0.92, 0.98),fontsize=fs) #(0.065, 0.98)
fig.legend(bbox_to_anchor=(.13, 1.08), loc='upper left', borderaxespad=0.,fontsize=10,framealpha=1)
#plt.title("Elevation of the 3 vortex from the 17 August to the 14 October",fontsize=fs)

#ax.set_facecolor('gainsboro')#whitesmoke gainsboro
plt.setp(ax.get_xticklabels(), rotation=45)
monthyearFmt = mdates.DateFormatter('%d %b')
ax.xaxis.set_major_formatter(monthyearFmt)
ax.tick_params(labelsize=fs)
ax2.tick_params(labelsize=fs)
ax2.text(-0.15,1.05,'c)',transform=ax2.transAxes,fontsize=16)
#_ = plt.xticks(rotation=45)

plt.savefig(savfile+'DATE_ALT.png',bbox_inches='tight',dpi=300,format='png')
#plt.show()