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
diroutput='/data/hlestrel/fig_bulle_'+initials+'/'
trackname='Ntrack-L1-'+initials
with gzip.open(diroutput+trackname,'rb') as f:
    infos = pickle.load(f)
# fenêtre de visualisation en classique
window_classic_lon=40
window_classic_lat=30
window_classic_alt=30

nametrackmodel = '/home/hlestrel/stageM1-vortex/track_model_CAN2017_BulleA_VO'
variable='VO'
mask_percent=1#0.25
suppr=False
arrangement=False
ano=False
lon_visu=30
lat_visu=25
alt_visu=20
# size of the lon/lat window
fenetre_lon=6 # fenêtre de recherche de la bulle
fenetre_lat=5
fenetre_alt=2
###### UPDATE TRACK
jour_butoir=27
print('jour butoir',jour_butoir)
print_irl=True
if os.path.isfile(nametrackmodel):
    update_track = True
    lon_update=True
    with gzip.open(nametrackmodel,'rb') as f:
        track_model = pickle.load(f)
    if suppr:
        print('suppr option engaged ! Deleting :',[idx, idx-1, idx-2])
        for xdel in range(257,264):
            del track_model[xdel]
    if arrangement: 
        idx=186 #À CHOISIR
    else:    
        idx=len(track_model)-1
    colori=(idx+200)*2
    ## EXP REVERSE #idx=113
    print('update_track = True ; continue track from the date ',track_model[idx]['date'].day,'/',track_model[idx]['date'].month,'/',track_model[idx]['date'].year)
    date00=track_model[idx]['date']
    if arrangement :  # À CHOISIR
        track_model[idx]['alt']=22
        track_model[idx]['lon']=10
        track_model[idx]['lat']=32
    altitude=track_model[idx]['alt']
    lon_mid=track_model[idx]['lon']
    if lon_mid<-180:
        lon_mid+=360
    lat_bulle=track_model[idx]['lat']
    if initials=="CAN2017":
        ## EXP REVERSE = -3
        ihours=3
    if initials=="AUS2009":
        ihours=3
else:
    update_track = False
    lon_update=False
    colori=600
    idx=0
    print('update_track = False ; new dico')
    dd0=infos[0]['date'] 
    altitude = infos[0]['mid'] 
    ihours=3
    track_model = {}
    if initials=="CAN2017":
        if nametrackmodel == '/home/hlestrel/stageM1-vortex/track_model_'+initials+'_BulleB2.2_test':
            dd0=datetime(2017,10,13,0)
            lon_mid=72
            lat_bulle =45
            iplus=0
            altitude = 21.4
        if nametrackmodel == '/home/hlestrel/stageM1-vortex/track_model_CAN2017_BulleA_VO':
            dd0=datetime(2017,8,29,3)
            lon_mid=40
            lat_bulle =54
            iplus=0
            altitude = 17.7
        else:
            ## Séparation B1 / B2
            dd0=datetime(2017,8,30,12)
            lon_mid=14
            lat_bulle =42
            iplus=0
            altitude = 18
        ## Séparation A / B
        # dd0=datetime(2017,8,24,12)
        # lon_mid=358
        # lat_bulle =54
        # iplus=0
        # altitude = 17
        # ##
    elif initials=="AUS2009":
        dd0=datetime(2009, 2, 10, 9)
        lon_mid=192
        lat_bulle=-46
        altitude=18.7
        iplus=0
        ihours=3
    else:
        lon_mid=infos[0]['lon_mid']
        lat_bulle =infos[0]['lat_mid']
        iplus=0
    isel0=0 
    date00=dd0
print(date00)
if (lon_mid+fenetre_lon/2)<359:
    loninf = lon_mid-fenetre_lon/2
    lonsup = lon_mid+fenetre_lon/2
elif (lon_mid+fenetre_lon/2)>359:
    loninf = lon_mid-fenetre_lon/2-359
    lonsup = lon_mid+fenetre_lon/2-359
if (lat_bulle-fenetre_lat/2)>-85 and (lat_bulle+fenetre_lat/2)<85:
    latinf = lat_bulle-fenetre_lat/2
    latsup = lat_bulle+fenetre_lat/2
elif (lat_bulle-fenetre_lat/2)<-85:
    latinf = -85
    latsup = lat_bulle+fenetre_lat/2
elif (lat_bulle+fenetre_lat/2)>85:
    latinf = lat_bulle-fenetre_lat/2
    latsup = 85
    
print('altitude :',altitude)
print('lon_mid :',lon_mid)
print('lat_mid :', lat_bulle)

dat = ECMWF('FULL-EA',date00,exp='VOZ')
if variable=='O3': 
    dat._get_var('O3')
elif variable=='T':
    dat._get_var('T')
elif variable=='VO':
    dat._get_var('VO')
else:
    print('NO VARIABLE DEFINED')
    exit()
# dat._get_var('U') 
# dat._get_var('V')
#dat._get_var('P')
#dat._get_var('PT')
try:
    dat._mkpscale() 
    dat._mkzscale()
    #dat._mkp() 
    #dat._mkthet() 
    #dat._mkpv()
    if update_track:
        if arrangement: 
            reponse=input('ARE YOU SURE ABOUT THE ARRANGEMENT OPTION ??? (y or n)')
            if reponse=='n':
                exit()
            print('ARRANGEMENT')
        else:
            ialts=track_model[idx]['ialts']
            ilons=track_model[idx]['ilons']
            ilats=track_model[idx]['ilats']
            min_O3=track_model[idx]['min_O3']
    else:
        iminalt=np.argmax(dat.attr['zscale']<altitude-fenetre_alt/2)
        ialt=np.argmax(dat.attr['zscale']<altitude)
        imaxalt=np.argmax(dat.attr['zscale']<altitude+fenetre_alt/2)
    O3mean = np.mean(dat.var[variable],axis=2)
    #Tmean = np.mean(dat.var['T'],axis=2)
    if ano:
        dat.var['O3ano'] = dat.var['O3'] - O3mean[:,:,np.newaxis]
        #dat.var['Tano'] = dat.var['T'] - Tmean[:,:,np.newaxis]
    if (loninf < 0):
        datp=dat.shift2west(150) 
        dats = datp.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
        if lat_bulle>0:
            dats00 = datp.extract(varss='All',latRange=(10,85))
        elif lat_bulle<0:
            dats00 = datp.extract(varss='All',latRange=(-80,-3))    
    else:
        dats = dat.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
        if lat_bulle>0:
            dats00 = dat.extract(varss='All',latRange=(10,85))
        elif lat_bulle<0:
            dats00 = dat.extract(varss='All',latRange=(-80,-3))

    if update_track==False:
        if ano:
            min_O3=np.amin(dats.var[variable+'ano'][imaxalt:iminalt])
        else: 
            min_O3=np.amin(dats.var[variable][imaxalt:iminalt])
        maskmax=min_O3+min_O3*mask_percent
        maskmin=min_O3-min_O3*mask_percent
        if ano:
            argwhere=np.argwhere(dats.var[variable+'ano'][imaxalt:iminalt]==min_O3)  
        else:
            argwhere=np.argwhere(dats.var[variable][imaxalt:iminalt]==min_O3)    
        ialts=imaxalt+argwhere[0][0]
        ilats=argwhere[0][1]
        ilons=argwhere[0][2]
    if arrangement:
        ialts=np.argmax(dats.attr['zscale']<altitude)
        track_model[idx]['ialts']=ialts
        ilons=np.argmax(dats.attr['lons']>lon_mid)-1
        track_model[idx]['ilons']=ilons
        ilats=np.argmax(dats.attr['lats']>lat_bulle)-1
        track_model[idx]['ilats']=ilats
        if ano:
            min_O3=dats.var[variable+'ano'][ialts][ilats][ilons]
        else:
            min_O3=dats.var[variable][ialts][ilats][ilons]
        track_model[idx]['min_O3']=min_O3
    print('min_var :',min_O3)
    
    # nilatmindats00 = np.argmax(85>latinf)-1
    # for ilat in range(len(dats00.var[variable][1])):
    #     for ilon in range(len(dats00.var[variable][1][1])):
    #         dats00.var[variable][ialts][ilat][ilon]=dats00.var[variable][ialts][ilat][ilon]-data_mean.var[variable][ialts][nilatmindats00+ilat]
    
    #dats.var[variable]=np.ma.masked_greater(dats.var[variable],maskmax)
    #dats.var[variable]=np.ma.masked_less(dats.var[variable],maskmin)
    
    if ano:
        ax=dats00.show(variable+'ano',ialts,cLines=False,projec=None,show=False,txt=str(variable)+' alt')
        ax2=dats00.chartlonz(variable+'ano',lat=dats.attr['lats'][ilats],levs=(ialts-20,ialts+20),show=False,txt=str(variable)+' lat')
    else:
        ax=dats00.show(variable,ialts,cLines=False,projec=None,show=False,txt=str(variable)+' alt')
        ax2=dats00.chartlonz(variable,lat=dats.attr['lats'][ilats],levs=(ialts-20,ialts+20),show=False,txt=str(variable)+' lat')
    jet= plt.get_cmap('binary') # jet
    colors = iter(jet(np.linspace(0,1,colori))) # attention à bien adapter au nombre de points que l'on souhaite tracer
    if update_track:
        for idxx in range(idx):
            altsi=track_model[idxx]['alt']
            lonsi=track_model[idxx]['lon']
            latsi=track_model[idxx]['lat']
            ax.scatter(lonsi+180,latsi,marker='o',linewidth=4,color=next(colors))
            if lonsi>0:
                ax2.scatter(lonsi,altsi,marker='o',linewidth=4,color=next(colors))
            else:
                ax2.scatter(lonsi+360,altsi,marker='o',linewidth=4,color=next(colors))
    else:
        if (lon_mid+fenetre_lon/2)<359:
            ax.scatter(dats.attr['lons'][ilons]-180,dats.attr['lats'][ilats],marker='o',linewidth=4,color=next(colors))
        elif (lon_mid+fenetre_lon/2)>359:
            ax.scatter(180+dats.attr['lons'][ilons],dats.attr['lats'][ilats],marker='o',linewidth=4,color=next(colors))
        ax2.scatter(dats.attr['lons'][ilons],dats.attr['zscale'][ialts],marker='o',linewidth=4,color=next(colors))
    fenetre_lat=2
    while date00.day !=jour_butoir: 
        idx+=1
    #for ihours in range(0,100,3):
        try:
            date00 = date00 + timedelta(hours=ihours)
            print("""
            ----------------------------------------
            """)
            print('date :',date00) 
            if lon_update:
                lon_update=False
            else:
                lon_mid=dats.attr['lons'][ilons]
                #if lon_mid-fenetre_lon/2>0 and (lon_mid+fenetre_lon/2)<359:
                if (lon_mid+fenetre_lon/2)<359:
                    loninf = lon_mid-fenetre_lon/2
                    lonsup = lon_mid+fenetre_lon/2
                # elif lon_mid-fenetre_lon/2<0:
                #     loninf = 360-lon_mid-fenetre_lon/2
                #     lonsup = lon_mid+fenetre_lon/2
                elif (lon_mid+fenetre_lon/2)>359:
                    loninf = lon_mid-fenetre_lon/2-359
                    lonsup = lon_mid+fenetre_lon/2-359
                
                lat_bulle =dats.attr['lats'][ilats]
                
                # if (lat_bulle-window_track_lat/2)>-85 or (lat_bulle+window_track_lat/2)<85:
                #     latinft = lat_bulle-window_track_lat/2
                #     latsupt = lat_bulle+window_track_lat/2
                # elif (lat_bulle-window_track_lat/2)<-85:
                #     latinft = -85
                #     latsupt = lat_bulle+window_track_lat/2
                # elif (lat_bulle+window_track_lat/2)>85:
                #     latinft = lat_bulle-window_track_lat/2
                #     latsupt = 85
                if (lat_bulle-fenetre_lat/2)>-85 and (lat_bulle+fenetre_lat/2)<85:
                    latinf = lat_bulle-fenetre_lat/2
                    latsup = lat_bulle+fenetre_lat/2
                elif (lat_bulle-fenetre_lat/2)<-85:
                    latinf = -85
                    latsup = lat_bulle+fenetre_lat/2
                elif (lat_bulle+fenetre_lat/2)>85:
                    latinf = lat_bulle-fenetre_lat/2
                    latsup = 85
                
                altitude = dats.attr['zscale'][ialts]
            print('altitude :',altitude)
            print('lon_mid :',lon_mid)
            print('lat_bulle :',lat_bulle)
            dat = ECMWF('FULL-EA',date00,exp='VOZ') 
            if variable=='O3': 
                dat._get_var('O3')
            elif variable=='T':
                dat._get_var('T')
            elif variable=='VO':
                dat._get_var('VO')
            else:
                print('NO VARIABLE DEFINED')
                exit()
            # dat._get_var('U')
            # dat._get_var('V')
            # dat._get_var('P')
            # dat._get_var('PT')
            try:
                dat._mkpscale() 
                dat._mkzscale()
                #dat._mkp() 
                #dat._mkthet() 
                #dat._mkpv()
                iminalt=np.argmax(dat.attr['zscale']<altitude-fenetre_alt/2)
                ialt=np.argmax(dat.attr['zscale']<altitude)
                imaxalt=np.argmax(dat.attr['zscale']<altitude+fenetre_alt/2)
                if ano:
                    O3mean = np.mean(dat.var[variable],axis=2)
                    dat.var[variable+'ano'] = dat.var[variable] - O3mean[:,:,np.newaxis]
                if (loninf < 0):
                    datp=dat.shift2west(150) 
                    dats = datp.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
                else:
                    dats = dat.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
                # data_mean=dat.zonal(vars='All')
                # #nilonmin = np.argmax(dat.attr['lons']>loninf)-1
                # #nilonmax = np.argmax(dat.attr['lons']>=lonsup)+1
                # nilatmin = np.argmax(dat.attr['lats']>latinf)-1
                # #nilatmax = np.argmax(dat.attr['lats']>latsup)
                # for altii in range(imaxalt,iminalt):
                #     for ilat in range(len(dats.var[variable][1])):
                #         for ilon in range(len(dats.var[variable][1][1])):
                #             dats.var[variable][altii][ilat][ilon]=dats.var[variable][altii][ilat][ilon]-data_mean.var[variable][altii][nilatmin+ilat]
                    
                maskmax=min_O3+min_O3*mask_percent
                maskmin=min_O3-min_O3*mask_percent
                if ano:
                    try:
                        if maskmin>maskmax:
                            mask1 = (dats.var[variable+'ano'][imaxalt:iminalt] < maskmin) & (dats.var[variable+'ano'][imaxalt:iminalt] > maskmax) 
                        elif maskmin<maskmax:
                            mask1 = (dats.var[variable+'ano'][imaxalt:iminalt] > maskmin) & (dats.var[variable+'ano'][imaxalt:iminalt] < maskmax) 
                        min_O3=np.amin(dats.var[variable+'ano'][imaxalt:iminalt][mask1])
                    except(TypeError,ValueError):
                        min_O3=np.amin(dats.var[variable+'ano'][imaxalt:iminalt])
                    print('min_var ;',min_O3)
                    argwhere=np.argwhere(dats.var[variable+'ano'][imaxalt:iminalt]==min_O3)  
                else:
                    try:
                        if maskmin>maskmax:
                            mask1 = (dats.var[variable][imaxalt:iminalt] < maskmin) & (dats.var[variable][imaxalt:iminalt] > maskmax) 
                        elif maskmin<maskmax:
                            mask1 = (dats.var[variable][imaxalt:iminalt] > maskmin) & (dats.var[variable][imaxalt:iminalt] < maskmax) 
                        min_O3=np.amin(dats.var[variable][imaxalt:iminalt][mask1])
                    except(TypeError,ValueError):
                        min_O3=np.amin(dats.var[variable][imaxalt:iminalt])
                    print('min_var ;',min_O3)
                    argwhere=np.argwhere(dats.var[variable][imaxalt:iminalt]==min_O3)
                ialts=imaxalt+argwhere[0][0]
                ilats=argwhere[0][1]
                ilons=argwhere[0][2]
                
                if date00.hour==999:
                    if (loninf < 0):
                        datsvisu= datp.extract(varss='All',lonRange=(dats.attr['lons'][ilons]-lon_visu/2,dats.attr['lons'][ilons]+lon_visu/2),latRange=(dats.attr['lats'][ilats]-lat_visu/2,dats.attr['lats'][ilats]+lat_visu/2))
                    else:
                        datsvisu= dat.extract(varss='All',lonRange=(dats.attr['lons'][ilons]-lon_visu/2,dats.attr['lons'][ilons]+lon_visu/2),latRange=(dats.attr['lats'][ilats]-lat_visu/2,dats.attr['lats'][ilats]+lat_visu/2))
                    # nilatmin = np.argmax(dat.attr['lats']>dats.attr['lats'][ilats]-lat_visu/2)-1
                    # for altii in range(ialts-6,ialts+6):
                    #     for ilat in range(len(datsvisu.var[variable][1])):
                    #         for ilon in range(len(datsvisu.var[variable][1][1])):
                    #             datsvisu.var[variable][altii][ilat][ilon]=datsvisu.var[variable][altii][ilat][ilon]-data_mean.var[variable][altii][nilatmin+ilat]
                    datsvisu.show(variable,ialts,cLines=False,projec=None,show=False)
                    datsvisu.show(variable,lat=dats.attr['lats'][ilats],levs=(ialts-6,ialts+6),show=False)
                ax2.scatter(dats.attr['lons'][ilons],dats.attr['zscale'][ialts],marker='o',linewidth=4,color=next(colors))
                if (dats.attr['lons'][ilons]+fenetre_lon/2)<359:
                    ax.scatter(dats.attr['lons'][ilons]-180,dats.attr['lats'][ilats],marker='o',linewidth=4,color=next(colors))
                elif (dats.attr['lons'][ilons]+fenetre_lon/2)>359:
                    ax.scatter(180+dats.attr['lons'][ilons],dats.attr['lats'][ilats],marker='o',linewidth=4,color=next(colors))
                gc.collect()
                track_model[idx] = {'date':date00,'alt':dats.attr['zscale'][ialts],'lon':dats.attr['lons'][ilons],'lat':dats.attr['lats'][ilats],'min_O3':min_O3,'ialts':ialts,'ilats':ilats,'ilons':ilons}
                with gzip.open(nametrackmodel,'wb') as f:
                    pickle.dump(track_model,f)
                if print_irl:
                    plt.pause(0.05)
            except('truc'):#(KeyError, ValueError):
                print('Donnée non lisible ou plsu gros problème.... ATTENTION !')
        except('problème'):#('IndexError'):
            print('IndexError')
            continue
    if print_irl==False:
        plt.show()
    with gzip.open(nametrackmodel,'wb') as f:
        pickle.dump(track_model,f)
except(KeyError):
    print('Première donnée non lisible')
    
    
### In case the suppr is needed :
# for xdel in range(625,len(track_model)):
#     del track_model[xdel]