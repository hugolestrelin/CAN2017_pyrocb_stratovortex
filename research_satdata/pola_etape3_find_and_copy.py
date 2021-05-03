#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Shows and store the data of the chosen sequences of CALIOP sections in function of the Cald dictionary, 'selCaliop_'+initials+'.pkl', from the directory defined in chemincald.
The chosen sequences had to be previously repertoried in a separate directory (diroutput).
It shows the sequences and allow the user to highlight the potential bubble/vortex with a rectangular seclector. This selector will allow to stored the useful data in order to use them in the meteorological model of the forth step of the finding_vortex program.
The user has just to run the program in order to select the bubble. When the figure is displayed, only the click is needed, and if the selection is completed, the '+' will store the data, and then display the next figure. If the user wants to stop it, he has to write 'e' on his keyboard.

input : chosen sequences (folderPath ; imagename) ; cald dictionary (chemincald)
output : track dictionary, in (diroutput)

@author: Bernard Legras, modified by Hugo Lestrelin

"""
from __future__ import print_function
from matplotlib.widgets import RectangleSelector
import numpy as np
import socket
if 'icare' in socket.gethostname():
    import sys,os
    sys.path.append('/home/hlestrelin/pylib')
from datetime import datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
from pyhdf.SD import SD, SDC
import glob
from astropy.time import Time, TimeDelta
from pyhdf import HDF, VS, V
import matplotlib.colors as colors
import geosat
from cartopy import feature
import cartopy.crs as ccrs
import argparse
from progress_bar import printProgressBar
import scipy.signal as ss

##
initials='CAN2017'
chemincald='/scratch/hlestrelin/fig_'+initials+'/selCaliop_'+initials+'.pkl' 
folderPath = '/scratch/hlestrelin/fig_'+initials+'/suivi_bulle/'
diroutput= folderPath
outfile = 'Ntrack-L1-'+initials
imagename='plot_'+initials+'.'

datasave=False 
iknown=None ## if first i known
Qs = 5.167*1.e-31
kbw = 1.0313

# find the data ; ICARE directories
if int(initials[-4:])>2019:
    dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v3.40'
    dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1.v3.40'
    dirCProf = '/DATA/LIENS/CALIOP/05kmCPro.v3.40'
else:
    dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v4.20'
    dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1.v4.20'
    dirCProf = '/DATA/LIENS/CALIOP/05kmCPro.v4.20'

##
# load dictionary of traces
with gzip.open(chemincald,'rb') as f:
    box,idate0,nbday,Cald = pickle.load(f)
    
# Check which traces are needed
cald_index=[]
for i in range(len(Cald)):
    fileplot=imagename+str(i)+'.png'
    nameplot = os.path.join( folderPath, fileplot )
    if os.path.isfile(nameplot) :
        cald_index.append(i)

# Does the update is needed or not :
track_index=[]
filetrack=outfile
nameplotcal = os.path.join( diroutput, filetrack )
if os.path.isfile(nameplotcal) :
    update = True
    with gzip.open(nameplotcal,'rb') as f:
        track = pickle.load(f)
    j=len(track)-1
    if track[j]['bot']==None:
        i=cald_index[j]
    else:
        i=cald_index[j+1]
    print('Update = True ; track read with ',len(track),' records ; starting from index '+str(i))

else:
    update = False
    track={}
    j=0
    i=cald_index[j]
    print('Update = False ; new track')

## Selectors : 
def line_select_callback(eclick, erelease):
    global dd,bot,mid,top,south,north,lon_mid,lat_mid
    'eclick and erelease are the press and release events'
    pos1 = np.where(lats<eclick.xdata)[0][0]
    tt = utc[pos1]
    dd = datetime.strptime(str(int(tt+20000000)),'%Y%m%d')+timedelta(days=tt-int(tt))
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    if eclick.ydata > erelease.ydata:
        top=eclick.ydata
        bot=erelease.ydata
    else:
        bot=eclick.ydata
        top=erelease.ydata
    mid=np.abs(eclick.ydata+erelease.ydata)/2
    if eclick.xdata < erelease.xdata:
        south=eclick.xdata
        north=erelease.xdata
    else:
        north=eclick.xdata
        south=erelease.xdata
    lat_mid=(eclick.xdata+erelease.xdata)/2
    pos_lon = np.where(lats<lat_mid)[0][0]
    lon_mid=lons[pos_lon]
    #track.append([i,Cald[i]['fname'],dd,bot,mid,top,south,north,lon_mid,lat_mid])
    print("Box defined : (lat : %3.2f, alt : %3.2f) --> (lat : %3.2f, alt : %3.2f)" % (x1, y1, x2, y2))
    print("Position of the cloud that will be saved for the day :",dd)
    print("alt_mid=",mid)
    print("lon_mid=",lon_mid)
    print("lat_mid=",lat_mid)

def toggle_selector(event):
    global breaker
    if event.key == '+':
        print('recording whole track and going to the next picture')
        plt.close('all')
    if event.key in ['e','E']:
        print('breaking and storing the results')
        breaker = True
        plt.close('all')
        

if iknown!=None:
    i=iknown

##
# Parameters of the plot
## CHANGEMENT TEMPORAIRE BOX
lat_min,lat_max,lon_min,lon_max=50,box[2],30,box[1] #°
#lat_min,lat_max,lon_min,lon_max=box[3],box[2],box[0],box[1] #°
year=Cald[1]['utc'].year
cm = 0
ext = [lon_min-180,lon_max-180,lat_min,lat_max]
xlims = [-180,180]
xlocs=None
itm = 20
yinf, ysup = 9, 24 
masked_alt_min=yinf
medfilter=False
widthfilter = 21
vmax = 20
breaker=False
point = True 
cmap = 'gist_ncar'

files={}
export={}
data_tot={}
printProgressBar(j, len(cald_index), prefix = 'Progrès de la selection', suffix = 'Complete', length = 50)
print('click and release in order to create a box around the bubble. You can then modify your box without consequences, until you press "+" add the data and pursue the selection or "e" to end it.')
while i <= cald_index[-1]:
    printProgressBar(j+1, len(cald_index), prefix = 'Progrès de la selection:', suffix = 'Complete', length = 50)
    date = datetime.strptime(Cald[i]['fname'][:10],'%Y-%m-%d')
    dirday = os.path.join(dirAProf,date.strftime('%Y/%Y_%m_%d'))
    dirdayL1 = os.path.join(dirL1,date.strftime('%Y/%Y_%m_%d'))
    dirdayC = os.path.join(dirCProf,date.strftime('%Y/%Y_%m_%d'))        
    try:
        if int(initials[-4:])>2019:
           file = os.path.join(dirday,'CAL_LID_L2_05kmAPro-Prov-V3-40.'+Cald[i]['fname']+'.hdf')
        else:
           file = os.path.join(dirday,'CAL_LID_L2_05kmAPro-Standard-V4-20.'+Cald[i]['fname']+'.hdf')
        hdf = SD(file,SDC.READ)
        hh = HDF.HDF(file,HDF.HC.READ)
        lons = hdf.select('Longitude').get()[Cald[i]['sel1'],1] % 360
        lats = hdf.select('Latitude').get()[Cald[i]['sel1'],1]
        meta = hh.vstart().attach('metadata')
        alts = np.array(meta.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
        a512 = np.ma.masked_less(hdf.select('Aerosol_Multiple_Scattering_Profile_532').get()[Cald[i]['sel1'],:],0)
        e512 = np.ma.masked_less(hdf.select('Extinction_Coefficient_532').get()[Cald[i]['sel1'],:],0)
        f512 = np.ma.masked_less(hdf.select('Aerosol_Layer_Fraction').get()[Cald[i]['sel1'],:],0)
        t512 = np.ma.masked_less(hdf.select('Total_Backscatter_Coefficient_532').get()[Cald[i]['sel1'],:],0)
        AProOK = True
        altii=0
        altiisup=0
        while alts[altii]>masked_alt_min:
            altii+=1
        while alts[altiisup]>ysup:
            altiisup+=1
    except:
        AProOK = False
        print('missing aerosol data')
    
    if point:
        fig, ax = plt.subplots()
        norm = colors.LogNorm(vmin=0.00001,vmax=0.01)
        im=ax.pcolormesh(lats,alts,t512.T,cmap=cmap,norm=norm) #[altiisup:altii]
        ax.set_ylim(yinf,ysup)
        ax.set_title('Total aerosol scattering ratio 512 ; index : '+str(i)+'/'+str(cald_index[-1]))
        plt.colorbar(im)
        dd=None; top = None; bot = None; mid = None; south = None; north = None ; lon_mid=None ; lat_mid=None
        toggle_selector.RS = RectangleSelector(ax, line_select_callback,
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True)
        cid=plt.connect('key_press_event', toggle_selector)
        plt.show()
        #track.append([i,Cald[i]['fname'],dd,bot,mid,top,south,north,lon_mid,lat_mid])
        #track[j]={'cald_index':i,'name':Cald[i]['fname'], 'date':dd, 'bot':bot, 'mid':mid, 'top':top, 'south':south, 'north':north, 'lon_mid':lon_mid, 'lat_mid':lat_mid}
        #export[j]={'lats':lats,'lons':lons,'alts':alts,'t512L1':t512L1_test,'malts':malts,'lbeta512_met':lbeta512_met,'mnd':mnd,'utc':utc}
        #files[j]={'hdf':hdf,'hh':hh}
        if breaker: break
          
    j+=1
    if j<len(cald_index):
        i=cald_index[j]
    else:
        break

if datasave:
    with gzip.open(diroutput+outfile,'wb') as f:
        pickle.dump(track,f)
    # with gzip.open(diroutput+'files','wb') as f:
    #     pickle.dump(files,f)
    with gzip.open(diroutput+'export','wb') as f:
        pickle.dump(export,f)
