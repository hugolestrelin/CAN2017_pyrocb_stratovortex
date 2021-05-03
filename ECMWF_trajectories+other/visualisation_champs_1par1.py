#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Montre et enregistre des sections du track choisi sur les autres champs que celui de l'O3 

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

initials='CAN2017'
diroutput='/data/hlestrel/fig_bulle_'+initials+'/visu_track_champs/'
nametrackmodel='/home/hlestrel/stageM1-vortex/track_model_CAN2017_BulleA_VO'

def toggle_selector(event):
    if event.key in ['c','C']:
        print('close all')
        plt.close('all')

inputest=False
loop=True
lon_visu=70
lat_visu=45
alt_visu=15
with gzip.open(nametrackmodel,'rb') as f:
    track_model = pickle.load(f)
print(track_model)
print('dernier indice = ',len(track_model))
idxi=np.int(input('indice du dictionnaire de départ  '))
for idx in range(idxi,len(track_model)-1,1):
    date00=track_model[idx]['date']
    min_O3=track_model[idx]['min_O3']
    if inputest:
        altitude=float(input('altitude ?  '))
        lon_mid=float(input('lon_mid ?  '))
        lat_bulle=float(input('lat_bulle ?  '))
    else:
        altitude=track_model[idx]['alt']
        lon_mid=track_model[idx]['lon']
        lat_bulle=track_model[idx]['lat']
        #ialts=track_model[idx]['ialts']
    if lon_mid<-180:
        lon_mid+=360
    lonrngmx=lon_mid+lon_visu/2
    lonrngmin=lon_mid-lon_visu/2
    if lonrngmx >= 359:
        lonrngmx+=-360
        lonrngmin+=-360
    print('idx', 'date00', 'altitude', 'lon_mid', 'lat_bulle')#,'min_O3')
    print(idx, date00, altitude, lon_mid, lat_bulle)#, min_O3)
    if os.path.isfile(diroutput+'dat_'+str(idx)) :
        with gzip.open(diroutput+'dat_'+str(idx),'rb') as f:
            datsvisu=pickle.load(f)
    else:
        dat = ECMWF('FULL-EA',date00,exp='VOZ') 
        dat._get_var('O3')
        dat._get_var('T')
        dat._get_var('U') 
        dat._get_var('V')
        dat._get_var('VO')
        dat._get_var('P')
        dat._get_var('PT')
        dat._get_var('Z')
        dat._mkpscale() 
        dat._mkzscale()
        dat._mkp() 
        dat._mkz()
        dat._mkthet() 
        dat._mkpv()
        O3mean = np.mean(dat.var['O3'],axis=2)
        Tmean = np.mean(dat.var['T'],axis=2)
        PVmean = np.mean(dat.var['PV'],axis=2)
        dat.var['PVano'] = dat.var['PV'] - PVmean[:,:,np.newaxis]
        dat.var['O3ano'] = dat.var['O3'] - O3mean[:,:,np.newaxis]
        dat.var['Tano'] = dat.var['T'] - Tmean[:,:,np.newaxis]
        if (lonrngmin < 0):
            datp=dat.shift2west(180) 
            datsvisu= datp.extract(varss='All',lonRange=(lonrngmin,lonrngmx),latRange=(lat_bulle-lat_visu/2,lat_bulle+lat_visu/2))
            #datsvisu= datp.extract(varss='All',lonRange=(lonrngmin,lonrngmx),latRange=(10,88))
        else:
            datsvisu = dat.extract(varss='All',lonRange=(lonrngmin,lonrngmx),latRange=(lat_bulle-lat_visu/2,lat_bulle+lat_visu/2))
        with gzip.open(diroutput+'dat_'+str(idx),'wb') as f:
            pickle.dump(datsvisu,f)
    # data_mean=dat.zonal(vars='All')
    # nilatmin = np.argmax(dat.attr['lats']>lat_bulle-lat_visu/2)-1
    # for ilat in range(len(datsvisu.var['O3'][1])):
    #     for ilon in range(len(datsvisu.var['O3'][1][1])): 
    #             datsvisu.var['O3'][ialts][ilat][ilon]=datsvisu.var['O3'][ialts][ilat][ilon]-data_mean.var['O3'][ialts][nilatmin+ilat]
    #             datsvisu.var['T'][ialts][ilat][ilon]=datsvisu.var['T'][ialts][ilat][ilon]-data_mean.var['T'][ialts][nilatmin+ilat]                             
    # ##
    # fenetre_alt=4
    # variable='O3ano'
    # mask_percent=0.55
    # iminalt=np.argmax(dat.attr['zscale']<altitude-fenetre_alt/2)
    # ialt=np.argmax(dat.attr['zscale']<altitude)
    # imaxalt=np.argmax(dat.attr['zscale']<altitude+fenetre_alt/2)
    # min_O3=np.amin(datsvisu.var[variable][imaxalt:iminalt])
    # maskmax=min_O3+min_O3*mask_percent
    # maskmin=min_O3-min_O3*mask_percent
    # if maskmin>maskmax:
    #     mask1 = (datsvisu.var[variable][imaxalt:iminalt] < maskmin) & (datsvisu.var[variable][imaxalt:iminalt] > maskmax) 
    # elif maskmin<maskmax:
    #     mask1 = (datsvisu.var[variable][imaxalt:iminalt] > maskmin) & (datsvisu.var[variable][imaxalt:iminalt] < maskmax) 
    # argwhere=np.argwhere(mask1==True)    
    # ialtsi=imaxalt+argwhere[:,0]
    # ilats=argwhere[:,1]
    # ilons=argwhere[:,2]
    # 
    # #conversion pour plot 3d selon projection georeferenced
    # R=6400
    # latitudes=(datsvisu.attr['lats'][ilats]-np.amin(datsvisu.attr['lats'][ilats]))*np.pi*R/180
    # longitudes=(datsvisu.attr['lons'][ilons]-np.amin(datsvisu.attr['lons'][ilons]))*np.pi*np.cos(np.deg2rad(datsvisu.attr['lats'][ilats]))*R/180
    # lon_mid_conv=(lon_mid-np.amin(datsvisu.attr['lons'][ilons]))*np.pi*np.cos(np.deg2rad(lat_bulle))*R/180
    # lat_bulle_conv=(lat_bulle-np.amin(datsvisu.attr['lats'][ilats]))*np.pi*R/180
    # aaa=np.meshgrid(datsvisu.attr['lons'],datsvisu.attr['lats']) 
    # #contour_lat_conv=(datsvisu.attr['lats']-np.amin(datsvisu.attr['lats'][ilats]))*np.pi*R/180
    # for i in range(len(aaa[1])):
    #     aaa[1][i]=(aaa[1][i]-np.amin(datsvisu.attr['lats'][ilats]))*np.pi*R/180
    # for i in range(len(datsvisu.attr['lats'])):
    #     aaa[0][i]=(datsvisu.attr['lons']-np.amin(datsvisu.attr['lons'][ilons]))*np.pi*np.cos(np.deg2rad(datsvisu.attr['lats'][i]))*R/180
    # lat_bulle_true=np.argwhere(datsvisu.attr['lats']==lat_bulle)   
    # contour_lon_conv=(datsvisu.attr['lons']-np.amin(datsvisu.attr['lons'][ilons]))*np.pi*np.cos(np.deg2rad(datsvisu.attr['lats'][lat_bulle_true]))*R/180     
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # img = ax.scatter(longitudes, latitudes, datsvisu.attr['zscale'][ialtsi], c=datsvisu.var[variable][imaxalt:iminalt][mask1], cmap='bone',alpha=0.1,linewidth=10,marker='*')
    # ax.scatter(lon_mid_conv,lat_bulle_conv,altitude,marker='H',linewidth=15,color='red')
    # ax.contourf(aaa[0],aaa[1],datsvisu.var[variable][ialt],cmap='plasma',offset=0) #ARG ZDIR NOT WORKING ???????!!!!
    # #bbb=datsvisu.var[variable][imaxalt-10:iminalt+10,lat_bulle_true, :]
    # #ccc=np.meshgrid(contour_lon_conv[0],datsvisu.attr['zscale_i'][imaxalt-10:iminalt+10])
    # #ax.set_zlim((10.,30.))
    # #ax2 = fig.add_subplot(122)
    # #cset=ax.contourf(ccc[0],ccc[1]*100,bbb[:,0,0,:],cmap='bone')
    # fig.colorbar(img)
    # plt.show()
    
    # fenetre_alt=4
    # variable='Tano'
    # iminalt=np.argmax(dat.attr['zscale']<altitude-fenetre_alt/2)
    # ialt=np.argmax(dat.attr['zscale']<altitude)
    # imaxalt=np.argmax(dat.attr['zscale']<altitude+fenetre_alt/2)
    # mask1 = datsvisu.var[variable][imaxalt:iminalt] < 999
    # argwhere=np.argwhere(mask1==True)    
    # ialts=imaxalt+argwhere[:,0]
    # ilats=argwhere[:,1]
    # ilons=argwhere[:,2]   
    # R=6400
    # latitudes=(datsvisu.attr['lats'][ilats]-np.amin(datsvisu.attr['lats'][ilats]))*np.pi*R/180
    # longitudes=(datsvisu.attr['lons'][ilons]-np.amin(datsvisu.attr['lons'][ilons]))*np.pi*np.cos(np.deg2rad(datsvisu.attr['lats'][ilats]))*R/180            
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # 
    # #img = ax.scatter(datsvisu.attr['lons'][ilons], datsvisu.attr['lats'][ilats], datsvisu.attr['zscale'][ialts], c=datsvisu.var[variable][imaxalt:iminalt], cmap=mymap,alpha=0.3)#plt.hot()
    # img = ax.scatter(longitudes, latitudes, datsvisu.attr['zscale'][ialts], c=datsvisu.var[variable][imaxalt:iminalt][mask1], cmap='plasma',alpha=0.4,linewidth=9,marker='*')
    # ax.scatter(lon_mid,lat_bulle,altitude,marker='H',linewidth=10,color='black')
    # fig.colorbar(img)
    # plt.connect('key_press_event', toggle_selector)
    # plt.show()
    dirout=diroutput
    ialts=np.argmax(datsvisu.attr['zscale']<altitude)
    if inputest:
        ialts=np.argmax(datsvisu.attr['zscale']<altitude)
    altmax=np.argmax(datsvisu.attr['zscale']<altitude+alt_visu*0.3)
    if altitude-alt_visu*0.7>0:
        altmin=np.argmax(datsvisu.attr['zscale']<altitude-alt_visu*0.7)
    else:
        altmin=np.argmax(datsvisu.attr['zscale']<=1)
    # # ax=datsvisu.show('PVano',ialts,projec=None,show=False,txt='Anomaly of potential vorticity at altitude '+str(round(altitude,2)))#,clim=(-0.1*10**-5,1.5*10**-5))#,savfile=dirout+'PValt-'+str(track_model[idx]['date'])+str(idx)+'.png'
    # # # if lon_mid>180:
    # # #     ax.scatter(lon_mid+180,lat_bulle,marker='+',linewidth=4,color='yellow')
    # # #     #ax3.scatter(lon_mid+180,lat_bulle,marker='+',linewidth=4,color='yellow')
    # # #     #ax5.scatter(lon_mid+180,lat_bulle,marker='+',linewidth=4,color='yellow')
    # # # else:
    # # 
    # # ax.scatter(lon_mid-180,lat_bulle,marker='+',linewidth=4,color='yellow')
    # # savfile=dirout+'PVanoalt-'+str(date00)+str(idx)+'.png'
    # # plt.savefig(savfile,bbox_inches='tight',dpi=300,format='png')
    # # ax2=datsvisu.chartlonz('PV',lat=lat_bulle,levs=(altmax,altmin),show=False,txt='Anomaly of potential vorticity')#,clim=(-0.1*10**-5,1.5*10**-5))#,savfile=dirout+'PVlat-'+str(track_model[idx]['date'])+str(idx)+'.png'
    # # if lon_mid>180:
    # #     ax2.scatter(lon_mid-360,altitude,marker='o',linewidth=4,color='yellow')
    # # else:
    # #     ax2.scatter(lon_mid,altitude,marker='o',linewidth=4,color='yellow')
    # # savfile2=dirout+'PVanolat-'+str(date00)+str(idx)+'.png'
    # # plt.savefig(savfile2,bbox_inches='tight',dpi=300,format='png')
    # # ax3=datsvisu.show('Tano',ialts,projec=None,show=False,txt='Tano alt '+str(altitude))#,savfile=dirout+'Tanoalt-'+str(track_model[idx]['date'])+str(idx)+'.png')#,clim=(-1,1))
    # # ax3.scatter(lon_mid-180,lat_bulle,marker='+',linewidth=4,color='yellow')
    # # savfile3=dirout+'Tanoalt-'+str(date00)+str(idx)+'.png'
    # # plt.savefig(savfile3,bbox_inches='tight',dpi=300,format='png')
    # # ax4=datsvisu.chartlonz('Tano',lat=lat_bulle,levs=(altmax,altmin),show=False,txt='Tano')#,savfile=dirout+'Tanolat-'+str(track_model[idx]['date'])+str(idx)+'.png')#,clim=(-1,1))
    # # if lon_mid>180:
    # #     ax4.scatter(lon_mid-360,altitude,marker='o',linewidth=4,color='yellow')
    # # else:
    # #     ax4.scatter(lon_mid,altitude,marker='o',linewidth=4,color='yellow')
    # # savfile4=dirout+'Tanolat-'+str(date00)+str(idx)+'.png'
    # # plt.savefig(savfile4,bbox_inches='tight',dpi=300,format='png')
    ax5=datsvisu.show('VO',ialts,scale=10**6,projec=None,show=False,txt=str('VO')+' alt '+str(altitude))#,savfile=dirout+'O3anoalt-'+str(track_model[idx]['date'])+str(idx)+'.png')#,clim=(-10*10**-7,-3*10**-7))
    ax5.scatter(lon_mid,lat_bulle,marker='+',linewidth=4,color='yellow')
    savfile5=dirout+'VOalt-'+str(date00)+str(idx)+'.png'
    plt.savefig(savfile5,bbox_inches='tight',dpi=300,format='png')
    ax6=datsvisu.chartlonz('VO',lat=lat_bulle,levs=(altmax,altmin),show=False,txt='VO')#,savfile=dirout+'O3anolat-'+str(track_model[idx]['date'])+str(idx)+'.png')#,clim=(-10*10**-7,-3*10**-7))
    if lon_mid>180:
        ax6.scatter(lon_mid-360,altitude,marker='o',linewidth=4,color='yellow')
    else:
        ax6.scatter(lon_mid,altitude,marker='o',linewidth=4,color='yellow')
    savfile6=dirout+'VOlat-'+str(date00)+str(idx)+'.png'
    plt.savefig(savfile6,bbox_inches='tight',dpi=300,format='png')
    plt.connect('key_press_event', toggle_selector)
    #plt.show()
    ## track_modelX[idx] = {'date':date00,'alt':altitude,'lon':lon_mid,'lat':lat_bulle,'min_O3':0,'ialts':0,'ilats':0,'ilons':0}
    ## with gzip.open(nametrackmodelX,'wb') as f:
    ##     pickle.dump(track_modelX,f)
    if loop==False:
        idx=np.int(input('indice du dictionnaire   '))
