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
nametrackmodel1 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_A.pkl'
nametrackmodel2 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_B1.pkl'
nametrackmodel3 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_B2.pkl'
savfile='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/'
track_model1 = pickle.load(open(nametrackmodel1,'rb'))
track_model2 = pickle.load(open(nametrackmodel2,'rb'))
track_model3 = pickle.load(open(nametrackmodel3,'rb'))
date00=track_model3['dates'][0]
with gzip.open(savfile+'dats00_PT_'+str(0),'rb') as f:
    dats00=pickle.load(f)
colori1=int(round(len(track_model1['alts'])/9)-2)
colori2=int(round(len(track_model2['alts'])/9)-2)
colori3=int(round(len(track_model3['alts'])/9)-2)

track_model1['PT']=[]
track_model2['PT']=[]
track_model3['PT']=[]
## chamngement Alt -> PT
for i in range(len(track_model1['dates'])):
    track_model1['PT'].append(dats00.var['PTmean'][np.where(dats00.attr['zscale']<=track_model1['alts'][i])[0][0],np.where(dats00.attr['lats']==np.int(np.round(track_model1['lats'][i])))[0][0],0])
for i in range(len(track_model2['dates'])):
    track_model2['PT'].append(dats00.var['PTmean'][np.where(dats00.attr['zscale']<=track_model2['alts'][i])[0][0],np.where(dats00.attr['lats']==np.int(np.round(track_model2['lats'][i])))[0][0],0])
for i in range(len(track_model3['dates'])):
    track_model3['PT'].append(dats00.var['PTmean'][np.where(dats00.attr['zscale']<=track_model3['alts'][i])[0][0],np.where(dats00.attr['lats']==np.int(np.round(track_model3['lats'][i])))[0][0],0])

jet1= plt.get_cmap('Reds')
colors12 = iter(jet1(np.linspace(0.15,1,colori1))) 
colors1 = iter(jet1(np.linspace(0.15,1,colori1))) 
colors12d = iter(jet1(np.linspace(0.1,0.15,10))) 
colors1d = iter(jet1(np.linspace(0.1,0.15,10))) 
jet2= plt.get_cmap('Blues')
colors22 = iter(jet2(np.linspace(0.15,1,colori2))) 
colors2 = iter(jet2(np.linspace(0.15,1,colori2)))
colors22d = iter(jet2(np.linspace(0.1,0.15,colori2))) 
colors2d = iter(jet2(np.linspace(0.1,0.15,colori2)))
jet3= plt.get_cmap('binary') #binary
colors32 = iter(jet3(np.linspace(0.15,1,colori3))) 
colors3 = iter(jet3(np.linspace(0.15,1,colori3)))
colors32d = iter(jet3(np.linspace(0.1,0.15,10))) 
colors3d = iter(jet3(np.linspace(0.1,0.15,10)))

fig, ax = plt.subplots(figsize=(11, 4))
plt.ylim(11,40)

ax2=ax.twinx()

lnwd=2
scission=30
ls='-'
mk1='o'
mk2='D'

## Pour la legende :
colorsX12 = jet1(np.linspace(0.15,0.85,2))
ax.plot(track_model1['dates'][0],track_model1['alts'][0],linewidth=lnwd,color=colorsX12[-1],label='Vortex A')#color=next(colors12)
colorsX22 = jet2(np.linspace(0.15,0.85,2))
ax.plot(track_model2['dates'][0],track_model2['alts'][0],linewidth=lnwd,color=colorsX22[-1],label='Vortex B1')#color=next(colors22)
colorsX32 = jet3(np.linspace(0.15,1,2))
ax.plot(track_model3['dates'][0],track_model3['alts'][0],linewidth=lnwd,color=colorsX32[-1],label='Vortex B2')#color=next(colors32)

itera=3
for i in range(itera,30,itera):
    
    ax.plot([track_model1['dates'][i-itera],track_model1['dates'][i]],[track_model1['alts'][i-itera],track_model1['alts'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=next(colors12d))
    ax2.plot([track_model1['dates'][i-itera],track_model1['dates'][i]],[track_model1['PT'][i-itera],track_model1['PT'][i]],linewidth=lnwd,marker=mk2,color=next(colors1d))#, linestyle='--')
    ax.plot([track_model2['dates'][i-itera],track_model2['dates'][i]],[track_model2['alts'][i-itera],track_model2['alts'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=next(colors22d))
    ax2.plot([track_model2['dates'][i-itera],track_model2['dates'][i]],[track_model2['PT'][i-itera],track_model2['PT'][i]],linewidth=lnwd,marker=mk2,color=next(colors2d))#, linestyle='--')
    ax.plot([track_model3['dates'][i-itera],track_model3['dates'][i]],[track_model3['alts'][i-itera],track_model3['alts'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=next(colors32d))
    ax2.plot([track_model3['dates'][i-itera],track_model3['dates'][i]],[track_model3['PT'][i-itera],track_model3['PT'][i]],linewidth=lnwd,marker=mk2,color=next(colors3d))#, linestyle='--')
  
  
itera=9
for i in range(scission+itera-3,len(track_model1['alts']),itera):
    ax.plot([track_model1['dates'][i-itera],track_model1['dates'][i]],[track_model1['alts'][i-itera],track_model1['alts'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=next(colors12))
    ax2.plot([track_model1['dates'][i-itera],track_model1['dates'][i]],[track_model1['PT'][i-itera],track_model1['PT'][i]],linewidth=lnwd,marker=mk2,color=next(colors1))#, linestyle='--')
for i in range(scission+itera-3,len(track_model2['alts']),itera):
    ax.plot([track_model2['dates'][i-itera],track_model2['dates'][i]],[track_model2['alts'][i-itera],track_model2['alts'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=next(colors22))
    ax2.plot([track_model2['dates'][i-itera],track_model2['dates'][i]],[track_model2['PT'][i-itera],track_model2['PT'][i]],linewidth=lnwd,marker=mk2,color=next(colors2))#, linestyle='--')
for i in range(scission+itera-3,len(track_model3['alts']),itera):
    ax.plot([track_model3['dates'][i-itera],track_model3['dates'][i]],[track_model3['alts'][i-itera],track_model3['alts'][i]],linewidth=lnwd, linestyle=ls,marker=mk1,color=next(colors32))
    ax2.plot([track_model3['dates'][i-itera],track_model3['dates'][i]],[track_model3['PT'][i-itera],track_model3['PT'][i]],linewidth=lnwd,marker=mk2,color=next(colors3))#, linestyle='--')
    #ax2.scatter(track_model3['lats'][i],track_model3['alts'][i],linewidth=lnwd,color=next(colors32))
    
ax.plot(track_model3['dates'][-1],track_model3['alts'][-1],linewidth=lnwd,color='k',marker='o',label='Baro altitude')#color=next(colors12)
ax2.plot(track_model3['dates'][-1],track_model3['PT'][-1],linewidth=lnwd,color='k',marker='D',label='Potential temperature')#color=next(colors22)
fs=18
ax.set(xlabel="Date",
       ylabel="Baro Altitude (km)") #baro
ax.xlabel_style = {'size': fs}
ax.ylabel_style = {'size': fs}
ax2.set_ylabel("Potential Temperature (K)") #baro
ax2.ylabel_style = {'size': fs}
#ax.legend(loc='upper left')
fig.legend(loc='upper left', bbox_to_anchor=(0.065, 0.98))
#plt.title("Elevation of the 3 vortex from the 17 August to the 14 October",fontsize=fs)

#ax.set_facecolor('gainsboro')#whitesmoke gainsboro
plt.setp(ax.get_xticklabels(), rotation=45)
monthyearFmt = mdates.DateFormatter('%d %b')
ax.xaxis.set_major_formatter(monthyearFmt)
#_ = plt.xticks(rotation=45)

plt.savefig(savfile+'DATE_ALT.png',bbox_inches='tight',dpi=300,format='png')
plt.show()