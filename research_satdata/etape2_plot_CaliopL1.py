#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Shows the sequence of CALIOP sections in function of the Cald dictionnary 'selCaliop_'+initials+'.pkl', from the directory defined in diroutput.

input : Cald directory
output : figures to be interpreted

After saving the figures, the ones that shows some king of bubble of dust have to be stored in a separated directory (like 'suivi_bulle') in the same directory as diroutput. This secondary directory will be needed to proceed the third step of the finding_vortex program. 

@author: Bernard Legras, modified by Hugo Lestrelin

"""
import pickle,gzip
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from datetime import datetime, timedelta
from pyhdf.SD import SD, SDC
from pyhdf import HDF, VS, V
import os
import numpy as np
import socket
if 'icare' in socket.gethostname():
    import sys,os
    sys.path.append('/home/hlestrelin/pylib')
import matplotlib.colors as colors
import geosat
import cartopy.crs as ccrs
import argparse
import scipy.signal as ss
from progress_bar import printProgressBar

initials = '++CAN2017'
diroutput = '/scratch/hlestrelin/fig_'+initials+'/'

if int(initials[-4:])>2019:
    dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v3.40'
    dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1.v3.40'
    dirCProf = '/DATA/LIENS/CALIOP/05kmCPro.v3.40'
else:
    dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v4.10'
    dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1.v4.10'
    dirCProf = '/DATA/LIENS/CALIOP/05kmCPro.v4.10'

i_known = None # To be changed if the first i is known
figsav = True
if figsav:
    print('Now the figures will be stored in the directory : '+diroutput)
show = False
vmax = 6 ## échelle ratio d'aérosols, le rendre adaptatif ?
medfilter = False
widthfilter = 81
option_antarctic = False # in case that the orbits are parralels to the longitudes (i.e around the poles)
cmap = 'gist_ncar'

Qs = 5.167*1.e-31
kbw = 1.0313

# load dictionary of traces
with gzip.open(diroutput+'selCaliop_'+initials+'.pkl','rb') as f:
    box,idate0,nbday,Cald = pickle.load(f)
    
# Does the update is needed or not :
cald_index=[]
for i in range(len(Cald)):
    fileplotcal='plot_'+initials+'.'+str(i)+'.png' # New file name defined here
    nameplotcal = os.path.join( diroutput, fileplotcal )
    if os.path.isfile(nameplotcal) :
        cald_index.append(i)
if len(cald_index)>1:
    update = True
    firsti = cald_index[-1]
else:
    update = False
    if i_known == None:
        firsti = 1
    else:
        firsti=i_known

lasti = len(Cald)
listi = range(firsti,lasti+1)

## Canada / Alaska
if initials == 'CAN2017' or initials == "CAN2011" or initials == '+CAN2017' or initials == '++CAN2017' :
    GEO = False
    cm = 0
    ext = [0-180,359-180,0,85]
    xlims = [-180,180]
    xlocs=None 
    itm = 20
    yinf = 9
    ysup = 30
## Australia
elif initials=='AUS2015' or initials=='AUS2014' or initials=='AUS2009':
    GEO = True
    cm = 0
    ext = [0-180,359-180,-85,0]
    xlims = [-180,180]
    xlocs = None
    itm=20
    yinf = 8
    ysup = 30

def toggle_selector(event):
	global breaker, cpr, i
	if event.key == 'z':
		print('recording picture in "/suivi_bulle"')
		#os.system('cp -r '+diroutput+name+' '+diroutput+'suivi_bulle')
		cpr=True
		plt.close('all')
	if event.key == 'd':
		print('next pic index '+str(i))
		plt.close('all')
	if event.key == 'q':
		print('previous pic, index ',i-1)
		plt.close('all')
		i=i-2
	if event.key == 's':
		print('breaking at '+str(i))
		plt.close('all')
		breaker = True
breaker=None
#for i in range(np.int(firsti),524):
#	name="plot_ALA2014."+str(i)+".png"
#	img = mpimg.imread(name)
#	plt.imshow(img)
#	cid=plt.connect('key_press_event', toggle_selector)
#	plt.show()
#	if breaker: break


## REGARDER CE QUI PREND DU TEMPS : COMPUTATION OU FIGURE ? ESSAYER RÉDUIRE COUT EN COUPANT À ALTITUDE
printProgressBar(0, lasti+1, prefix = "Progrès de l'affichage", suffix = 'Complete', length = 50)
for i in listi:
    cpr=None
    # generate date and daily directory name
    date = datetime.strptime(Cald[i]['fname'][:10],'%Y-%m-%d')
    dirday = os.path.join(dirAProf,date.strftime('%Y/%Y_%m_%d'))
    dirdayL1 = os.path.join(dirL1,date.strftime('%Y/%Y_%m_%d'))
    dirdayC = os.path.join(dirCProf,date.strftime('%Y/%Y_%m_%d'))

    print('i',i)
    print(Cald[i]['fname'])
    #try:
        #if int(initials[-4:])>2019:
        #    file = os.path.join(dirday,'CAL_LID_L2_05kmAPro-Prov-V3-40.'+Cald[i]['fname']+'.hdf')
        #else:
        #    file = os.path.join(dirday,'CAL_LID_L2_05kmAPro-Standard-V4-10.'+Cald[i]['fname']+'.hdf')
        #hdf = SD(file,SDC.READ)
        #hh = HDF.HDF(file,HDF.HC.READ)
        #lons = hdf.select('Longitude').get()[Cald[i]['sel1'],1] % 360
        #lats = hdf.select('Latitude').get()[Cald[i]['sel1'],1]
        #meta = hh.vstart().attach('metadata')
        #alts = np.array(meta.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
        #a512 = np.ma.masked_less(hdf.select('Aerosol_Multiple_Scattering_Profile_532').get()[Cald[i]['sel1'],:],0)
        #e512 = np.ma.masked_less(hdf.select('Extinction_Coefficient_532').get()[Cald[i]['sel1'],:],0)
        #f512 = np.ma.masked_less(hdf.select('Aerosol_Layer_Fraction').get()[Cald[i]['sel1'],:],0)
        #t512 = np.ma.masked_less(hdf.select('Total_Backscatter_Coefficient_532').get()[Cald[i]['sel1'],:],0)
        #AProOK = True
    #except:
        #AProOK = False
        #print('missing aerosol data')

    #try:
        #if int(initials[-4:])>2019:
        #    fileC = os.path.join(dirdayC,'CAL_LID_L2_05kmCPro-Prov-V3-40.'+Cald[i]['fname']+'.hdf')
        #else:
        #    fileC = os.path.join(dirdayC,'CAL_LID_L2_05kmAPro-Standard-V4-10.'+Cald[i]['fname']+'.hdf')
        #hdfC = SD(fileC,SDC.READ)
        #hhC = HDF.HDF(fileC,HDF.HC.READ)
        #t512C = np.ma.masked_less(hdfC.select('Total_Backscatter_Coefficient_532').get()[Cald[i]['sel1'],:],0)
        #metaC = hhC.vstart().attach('metadata')
        #altsC = np.array(metaC.read()[0][meta.field('Lidar_Data_Altitudes')._idx])
        #CProOK = True
    #except:
        #CProOK = False
        #print('missing cloud data')

    try:
        if int(initials[-4:])>2019:
            fileL1 = os.path.join(dirdayL1,'CAL_LID_L1-ValStage1-V3-40.'+Cald[i]['fname']+'.hdf')
        else:
            fileL1 = os.path.join(dirdayL1,'CAL_LID_L1-Standard-V4-10.'+Cald[i]['fname']+'.hdf')
        hdf1 = SD(fileL1,SDC.READ)
        hh1 = HDF.HDF(fileL1,HDF.HC.READ)
        sel1L1 = Cald[i]['sel1L1'][:,0]

        t512L1 = hdf1.select('Total_Attenuated_Backscatter_532').get()[sel1L1,:]
        t512L1 = np.ma.masked_less(t512L1,0)
        lons1 = hdf1.select('Longitude').get()[sel1L1].flatten() % 360
        lats1 = hdf1.select('Latitude').get()[sel1L1].flatten()
        mnd1 = hdf1.select('Molecular_Number_Density').get()[sel1L1,:]
        lbeta512_met = np.log(1000 * mnd1 * Qs / (kbw*8*np.pi/3))
        meta1 = hh1.vstart().attach('metadata')
        alts1 = np.array(meta1.read()[0][meta1.field('Lidar_Data_Altitudes')._idx])
        meta1 = hh1.vstart().attach('metadata')
        malts1 = np.array(meta1.read()[0][meta1.field('Met_Data_Altitudes')._idx])
        # calculation of the molecular backscatter
        lbeta512_lid = np.empty(shape=t512L1.shape)
        for jy in range(len(lats1)):
            lbeta512_lid[jy,:] = np.interp(alts1,malts1[::-1],lbeta512_met[jy,::-1])
        if medfilter:
            sr512raw = t512L1/np.exp(lbeta512_lid)
            sr512= ss.medfilt(sr512raw,kernel_size=(widthfilter,1))
        else:
            sr512 = t512L1/np.exp(lbeta512_lid)
        L1OK = True
    except:
        print ('no L1 data')
        print (fileL1)
        L1OK = False

#%% Here we get the himawari scene which is the closest for this time

    GEOOK = False
    if GEO:
        try:
            utc = Cald[i]['utc']
            utc -= timedelta(minutes=(utc.minute%itm))
            try:
                ah = satload(utc)
            except:
                try:
                    ah = satload(utc+timedelta(minutes=itm))
                except:
                    try:
                        ah = satload(utc+timedelta(minutes=2*itm))
                    except:
                        ah = satload(utc-timedelta(minutes=itm))
            ah._get_IR0()
            ph = geosat.SatGrid(ah,gg)
            ph._sat_togrid('IR0')
            GEOOK = True
        except:
            GEOOK = False
    fig = plt.figure(figsize=(6,8))
    if L1OK:
        plt.subplot(2,1,1)
        norm1 = colors.LogNorm(vmin=0.00001,vmax=0.1)
        print('dims',lats1.shape,alts1.shape,sr512.shape)
        if option_antarctic:
            plt.pcolormesh(lons1,alts1,sr512.T,cmap=cmap,vmin=0,vmax=6)
        else:
            plt.pcolormesh(lats1,alts1,sr512.T,cmap=cmap,vmin=0,vmax=6)
        plt.title('{:d} {:16} {:.2f} {:.2f}'.format(i,Cald[i]['utc'].strftime('%Y-%m-%dT%H:%M'),lons1[0],lons1[-1]))
        plt.ylim(yinf,ysup)
        plt.title('Scattering ratio 512')
        plt.colorbar()
        #cid=plt.connect('key_press_event', toggle_selector)

    # if AProOK:
    #    plt.subplot(2,2,2)
    #    norm = colors.LogNorm(vmin=0.00001,vmax=0.01)
    #    plt.pcolormesh(lats,alts,t512.T,cmap=cmap,norm=norm)
    #    plt.ylim(yinf,ysup)
    #    plt.title('Total scattering aerosols 512')
    #    plt.colorbar()
    # if CProOK:
    #    plt.subplot(2,2,3)
    #    norm = colors.LogNorm(vmin=0.001,vmax=0.1)
    #    plt.pcolormesh(lats,altsC,t512C.T,cmap=cmap,norm=norm)
    #    plt.ylim(yinf,ysup)
    #    plt.title('Total scattering clouds 512')
    #    plt.colorbar()
    tran = ccrs.PlateCarree(central_longitude=cm)
    proj = tran
    # geostationary plot
    ax = plt.subplot(2,1,2,projection = proj)
    ax.set_extent(box,ccrs.PlateCarree())
    if GEO & GEOOK:
        buf = ph.var['IR0'][subgg.corner[1]:subgg.corner[1]+subgg.box_biny,
                            subgg.corner[0]:subgg.corner[0]+subgg.box_binx]
        ax.imshow(buf,extent=ext,cmap='jet',clim=(190,300),
              transform=tran,interpolation='nearest',origin='lower')
    ax.coastlines('50m')
    lons1[lons1>180] -= 360
    # to make sure, plot it twice
    ax.plot(lons1+cm,lats1,'k')
    ax.plot(lons1-cm,lats1,'k')
    ax.set_xlim(xlims[0],xlims[1])
    gl = ax.gridlines(draw_labels=True,xlocs=xlocs,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    if GEO & GEOOK: plt.title(utc.strftime('Hima %Y-%m-%d %H:%M'))
    fig.suptitle(Cald[i]['fname'])
    if figsav:
        plt.savefig(diroutput+'plot_'+initials+'.'+str(i)+'.png',dpi=300,bbox_inches='tight')
        print('The figures has been stored in the directory : '+diroutput)
        #if cpr:
        #        import os 
        #        os.system('cp -r '+diroutput+'plot_'+initials+'.'+str(i)+'.png '+diroutput+'suivi_bulle')
    if show: plt.show()
    plt.close(fig=fig)
    #if AProOK:
       #hh.close()
       #hdf.end()
    #if CProOK:
       #hhC.close()
       #hdfC.end()
    if L1OK:
        hh1.close()
        hdf1.end()
    printProgressBar(i, lasti+1, prefix = "Progrès de l'affichage", suffix = 'Complete', length = 50)
    if breaker: break

## REGARDER POSSIBILITÉ COMMANDE UNIX SUR PYTHON
