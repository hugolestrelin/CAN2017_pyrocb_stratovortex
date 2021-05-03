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
savfile='/data/hlestrel/fig_bulle_'+initials+'/fig_publi/VA_portion/min_O3/'
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
altitude=track_model3['z'][idebut]
lon_mid=track_model3['lons'][idebut]
lat_bulle=track_model3['lats'][idebut]
jour_butoir=26
date00=track_model3['dates'][idebut] #choix de l'indice pour l'image de fond
ihours=3
i=0
while date00.day !=jour_butoir: 
    date00 = date00 + timedelta(hours=ihours)
    print(date00,i)

    namedats=savfile+'dats_'+str(i)
    if os.path.isfile(namedats):
        with gzip.open(namedats,'rb') as f:
            dats00=pickle.load(f)
    else:    
        dat = ECMWF('FULL-EA',date00,exp='VOZ') 
        dat._get_var('O3')
        dat._mkpscale() 
        dat._mkzscale()
        #dat._mkp() 
        #dat._mkthet() 
        #dat._mkpv()
        O3mean = np.mean(dat.var['O3'],axis=2)
        dat.var['O3ano'] = dat.var['O3'] - O3mean[:,:,np.newaxis]
        datp=dat.shift2west(140) 
        dats00 = datp.extract(varss=['O3','O3ano'],lonRange=(lonrangemin,lonrangemax),latRange=(latrangemin,latrangemax))
            
        ## Enregistrement des données à réutiliser pour ce plot-ci : 
        with gzip.open(namedats,'wb') as f:
            pickle.dump(dats00,f)
    i+=1