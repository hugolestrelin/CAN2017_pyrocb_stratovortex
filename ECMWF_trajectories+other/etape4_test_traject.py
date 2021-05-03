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

#### AJOUTER COMPOSITE ET TRACE DE LA SELCTION (e.g point en mid)
classic=False 
figsave=False 
show=True
mean_PV=False
figsavemean=False 
showmean=False
composite=False
trajectoire=False
visualization_track=False
print_traj=True
initials='CAN2017'
diroutput='/data/hlestrel/fig_bulle_'+initials+'/'
trackname='Ntrack-L1-'+initials
with gzip.open(diroutput+trackname,'rb') as f:
    infos = pickle.load(f)
# fenêtre de visualisation en classique
window_classic_lon=40
window_classic_lat=30
window_classic_alt=30

# Does the update is needed or not :
update_index=[]
sel3=list(range(0,25,3))
filetrack='plotVPalt-'+str(infos[0]['cald_index'])+'.png'
nameplotcal = os.path.join( diroutput, filetrack )
if os.path.isfile(nameplotcal) :
    update = True
    for i in range(len(infos)):
        filetrack='plotVPalt-'+str(infos[i]['cald_index'])+'.png'
        nameplotcal = os.path.join( diroutput, filetrack )
        if os.path.isfile(nameplotcal) :
            update_index.append(i)
            firsti=update_index[-1]+1
    print('Update = True ; continue plots from index '+str(infos[firsti]['cald_index'])+' / '+str(infos[len(infos)-1]['cald_index']))
else:
    update = False
    firsti=0
    print('Update = False ; new plots')

if classic:
    printProgressBar(firsti, len(infos), prefix = 'Progrès du plot', suffix = 'Complete', length = 50)
    for i in range(firsti,len(infos)):
        printProgressBar(i, len(infos), prefix = 'Progrès du plot', suffix = 'Complete', length = 50)
        dd=infos[i]['date']
        sel3=list(range(0,25,3))
        isel=0
        while sel3[isel]<dd.hour:
            isel+=1
        if (dd.hour!=sel3[isel]):
            if dd.hour<21:
                diff=np.abs(sel3[isel]-dd.hour)
                hour=dd.hour+diff
                date = datetime(dd.year, dd.month, dd.day, hour)
            else:
                hour=sel3[isel-1]
                date = datetime(dd.year, dd.month, dd.day, hour)
        else:
            date=datetime(dd.year, dd.month, dd.day, dd.hour)
        
        
        lon_mid=infos[i]['lon_mid']
        #if lon_mid<(180-fenetre_lon/2):
        loninf = lon_mid-window_classic_lon/2
        lonsup = lon_mid+window_classic_lon/2
        #else:
            #loninf = (lon_mid-(360-fenetre_lon/2))-fenetre_lon/2
            #lonsup = (lon_mid-(360-fenetre_lon/2))+fenetre_lon/2
            
        lat_bulle =infos[i]['lat_mid']
        if (lat_bulle-window_classic_lat/2)>-85 or (lat_bulle+window_classic_lat/2)<85:
            latinf = lat_bulle-window_classic_lat/2
            latsup = lat_bulle+window_classic_lat/2
        elif (lat_bulle-window_classic_lat/2)<-85:
            latinf = -85
            latsup = lat_bulle+window_classic_lat/2
        elif (lat_bulle+window_classic_lat/2)>85:
            latinf = lat_bulle-window_classic_lat/2
            latsup = 85
        
        altitude = infos[i]['mid']
        
        dat = ECMWF('FULL-EA',date,exp='VOZ') 
        #dat._get_var('U') 
        #dat._get_var('V')
        dat._get_var('T')
        dat._get_var('O3') 
        #dat._get_var('P')
        #dat._get_var('PT') 
        try:
            dat._mkpscale() 
            dat._mkzscale()
            #dat._mkp() 
            #dat._mkthet() 
            #dat._mkpv()
            ialt=0
            while dat.attr['zscale'][ialt]>altitude:
                ialt+=1
            if (loninf < 0):
                datp=dat.shift2west(150) 
                dats = datp.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
            else:
                dats = dat.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
            
            data_mean=dat.zonal(vars='All')
            
            for altii in range(len(dats.var['O3'])):
                for ilat in range(len(dats.var['O3'][1])):
                    for ilon in range(len(dats.var['O3'][1][1])):
                        dats.var['O3'][altii][ilat][ilon]=dats.var['O3'][altii][ilat][ilon]-data_mean.var['O3'][altii][ilat]
                        dats.var['T'][altii][ilat][ilon]=dats.var['T'][altii][ilat][ilon]-data_mean.var['T'][altii][ilat]
                        
            if figsave:
                dats.show('O3',ialt,cLines=False,projec='nearside',show=False,savfile=diroutput+'plotO3alt-'+str(infos[i]['cald_index'])+'.png',txt='O3 alt = '+str(altitude)+'km '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
                dats.show('T',ialt,cLines=False,projec='nearside',show=False,savfile=diroutput+'plotTalt-'+str(infos[i]['cald_index'])+'.png',txt='Temperature anomaly alt = '+str(altitude)+'km '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
                plt.close('all')
                dats.chartlonz('O3',lat=lat_bulle,levs=(ialt-10,ialt+10),show=False,savfile=diroutput+'plotO3lat-'+str(infos[i]['cald_index'])+'.png',txt='O3'+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
                plt.close('all')
                dats.chartlonz('T',lat=lat_bulle,levs=(ialt-10,ialt+10),show=False,savfile=diroutput+'plotTlat-'+str(infos[i]['cald_index'])+'.png',txt='Temperature anomaly '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
                plt.close('all')
                print('The plot for the index '+str(infos[i]['cald_index'])+' has been saved')
            if show:
                dats.show('O3',ialt,cLines=False,projec='nearside',txt='Ano O3 @ alt = %3.2f km '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour) % (altitude))
                dats.show('T',ialt,cLines=False,projec='nearside',txt='Ano T @ alt = %3.2f km '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour) % (altitude))
                dats.chartlonz('O3',lat=lat_bulle,levs=(ialt-10,ialt+10),txt='Ano O3 @ lat = %3.2f km'+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour) % (lat_bulle))
                dats.chartlonz('T',lat=lat_bulle,levs=(ialt-10,ialt+10),txt='Temperature ano @ lat = %3.2f'+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour) % (lat_bulle))
        
            if mean_PV:
                # Test pour récupérer et déduire la moyenne de la PV !!! Problème de mémoire, pas possible.
                # FAIRE UN TRUC POUR ENREGISTRER DONNÉES LORS 1ERE FOIS PUIS UPDATE...
                nbday=22 #nbr de jour sur lesquels s'effectue la moeynne, ne pas changer ! 
                date0= date - timedelta(days=int(nbday/2))
                # if date0.month==datefirst.month:
                #     if date0.day<=datefirst.day:
                #         date0=datefirst
                # elif date0.month<datefirst.month:
                #     date0=datefirst
                datsitrack='datsi_'+str(date0.year)+'-'+str(date0.month)+'-'+str(date0.day)+'_'+str(loninf)+'_'+str(lonsup)+'_'+str(latinf)+'_'+str(latsup)
                namedatsi = os.path.join( diroutput, datsitrack )
                if os.path.isfile(namedatsi) :
                    print('Date and place already known for the mean')
                    with gzip.open(diroutput+datsitrack,'rb') as f:
                        datsi = pickle.load(f)
                else:
                    print('Date or place not known for the mean, compute and save...')
                    dati=ECMWF('FULL-EA',date0,exp='VOZ')
                    dati._get_var('U') 
                    dati._get_var('V')
                    dati._get_var('T')
                    dati._get_var('VO') 
                    dati._get_var('P')
                    dati._get_var('PT') 
                    #try:
                    dati._mkpscale() 
                    dati._mkzscale()
                    dati._mkp() 
                    dati._mkthet() 
                    dati._mkpv()
                    if (loninf < 0):
                        datpi=dati.shift2west(150) 
                        datsi = datpi.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
                    else:
                        datsi = dati.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
                    with gzip.open(diroutput+datsitrack,'wb') as f:
                        pickle.dump(datsi,f)
                    #except:
                        #print('Problème avec les données du jour :',date)
                printProgressBar(1, nbday-1, prefix = 'Progrès du calcul de la moyenne', suffix = 'Complete', length = 50)
                for iday in range(nbday-1):
                    datei = date0 + timedelta(days=iday+1)
                    print(datei)
                    print(i)
                    datsitrack='datsi_'+str(datei.year)+'-'+str(datei.month)+'-'+str(datei.day)+'_'+str(loninf)+'_'+str(lonsup)+'_'+str(latinf)+'_'+str(latsup)
                    namedatsi = os.path.join( diroutput, datsitrack )
                    if os.path.isfile(namedatsi) :
                        print('Date and place already known for the mean')
                        with gzip.open(diroutput+datsitrack,'rb') as f:
                            datsi = np.append(datsi,pickle.load(f))
                    else:
                        print('Date or place not known for the mean, compute and save...')
                        dati=ECMWF('FULL-EA',datei,exp='VOZ')
                        dati._get_var('U') 
                        dati._get_var('V')
                        dati._get_var('T')
                        dati._get_var('VO') 
                        dati._get_var('P')
                        dati._get_var('PT') 
                        try:
                            dati._mkpscale() 
                            dati._mkzscale()
                            dati._mkp() 
                            dati._mkthet() 
                            dati._mkpv()
                            if (loninf < 0):
                                datpi=dati.shift2west(150) 
                                datsi = np.append(datsi,datpi.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)))
                            else:
                                datsi = np.append(datsi,dati.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup))) 
                            with gzip.open(diroutput+datsitrack,'wb') as f:
                                pickle.dump(datsi,f)
                            printProgressBar(iday, nbday-1, prefix = 'Progrès du calcul de la moyenne', suffix = 'Complete', length = 50)
                        except('U not found or read error'):
                            print('Problème avec les données du jour :',datei)
                            continue
                        
                for i in range(len(datsi)-1):
                    datsi[0].var['PV']+=datsi[i+1].var['PV']
                
                datsi[0].var['PV']=(datsi[0].var['PV']/len(datsi))*10**6
                datsi[0].var['PV']=dats.var['PV']-datsi[0].var['PV']
            
                if figsavemean:
                    datsi[0].show('PV',ialt,cLines=False,projec='nearside',show=False,savfile=diroutput+'plotVPaltmean-'+str(infos[i]['cald_index'])+'.png',txt='Anomaly of the potential vorticity alt = '+str(altitude)+'km '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
                    plt.close('all')
                    datsi[0].chartlonz('PV',lat=lat_bulle,levs=(ialt-10,ialt+10),show=False,savfile=diroutput+'plotVPlatmean-'+str(infos[i]['cald_index'])+'.png',txt='Anomaly of the potential vorticity '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
                    plt.close('all')
                    print('The mean plot for the index '+str(infos[i]['cald_index'])+' has been saved')
                if showmean:
                    datsi[0].show('PV',ialt,cLines=False,projec='nearside',txt='Anomaly of the potential vorticity alt = %3.2f km '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour) % (altitude))
                    datsi[0].chartlonz('PV',lat=lat_bulle,levs=(ialt-10,ialt+10),txt='Anomaly of the potential vorticity '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
                    
        except(KeyError):
            print('Problème avec les données de l"index :',infos[i]['cald_index'])
            continue
        
        
#### AJOUTER LA FONCTION COMPOSITE et trajectoire :
if composite:
    a=1

nametrackmodel = '/home/hlestrel/stageM1-vortex/track_model_'+initials+'_BulleB2.2_test'
#### Si point ne fonctionne pas, rechercher "carré" de minimum avce point en son barycentre
if trajectoire:
    nametrackmodel = '/home/hlestrel/stageM1-vortex/track_model_'+initials+'_BulleB2.2_test'
    variable='O3'
    mask_percent=0.25
    suppr=False
    lon_visu=30
    lat_visu=25
    alt_visu=20
    # size of the lon/lat window
    fenetre_lon=6 # fenêtre de recherche de la bulle
    fenetre_lat=5
    fenetre_alt=2
    ###### UPDATE TRACK
    jour_butoir=30
    colori=1800
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
        ### arrangement : 
        #idx=à choisir
        idx=len(track_model)-3
        ## EXP REVERSE #idx=113
        print('update_track = True ; continue track from the date ',track_model[idx]['date'].day,'/',track_model[idx]['date'].month,'/',track_model[idx]['date'].year)
        date00=track_model[idx]['date']
        #if jour_butoir in track_model[:]['date'].day: sys.exit('jour_butoir deja fait')
        ## Arrangement : 
        #track_model[idx]['alt']=21.4
        #track_model[idx]['lon']=25
        #track_model[idx]['lat']=48
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
        while sel3[isel0]<dd0.hour: 
            isel0+=1 
        if (dd0.hour!=sel3[isel0]):
            if dd0.hour<21: 
                diff0=np.abs(sel3[isel0]-dd0.hour)
                hour0=dd0.hour+diff0
                date00 = datetime(dd0.year, dd0.month, dd0.day+iplus, hour0)
            else:
                hour0=sel3[isel0-1]
                date00 = datetime(dd0.year, dd0.month, dd0.day+iplus, hour0)
        else:
            date00=datetime(dd0.year, dd0.month, dd0.day+iplus, dd0.hour)
    print(date00)
    # loninft = lon_mid-window_track_lon/2
    # lonsupt = lon_mid+window_track_lon/2
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
       
    print('altitude :',altitude)
    print('lon_mid :',lon_mid)
    print('lat_mid :', lat_bulle)
    
    dat = ECMWF('FULL-EA',date00,exp='VOZ') 
    dat._get_var(variable)
    #dat._get_var('T')
    #dat._get_var('U') 
    #dat._get_var('V')
    #dat._get_var('VO')
    #dat._get_var('P')
    #dat._get_var('PT')
    try:
        dat._mkpscale() 
        dat._mkzscale()
        #dat._mkp() 
        #dat._mkthet() 
        #dat._mkpv()
        if update_track:
            ## if Arrangement : à commenter
            #print('truc')
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
        dat.var[variable+'ano'] = dat.var[variable] - O3mean[:,:,np.newaxis]
        #dat.var['Tano'] = dat.var['T'] - Tmean[:,:,np.newaxis]
        if (loninf < 0):
            datp=dat.shift2west(150) 
            dats = datp.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
            if lat_bulle>0:
                dats00 = datp.extract(varss='All',latRange=(10,85))
            elif lat_bulle<0:
                dats00 = datp.extract(varss='All',latRange=(-80,-3))
            #datsvisu= datp.extract(varss='All',lonRange=((lon_mid-lon_visu/2),(lon_mid+lon_visu/2)),latRange=(lat_bulle-lat_visu/2,lat_bulle+lat_visu/2))
        else:
            dats = dat.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
            if lat_bulle>0:
                dats00 = dat.extract(varss='All',latRange=(10,85))
            elif lat_bulle<0:
                dats00 = dat.extract(varss='All',latRange=(-80,-3))
            #datsvisu = dat.extract(varss='All',lonRange=(lon_mid-lon_visu/2,lon_mid+lon_visu/2),latRange=(lat_bulle-lat_visu/2,lat_bulle+lat_visu/2))
    
        
        # data_mean=dat.zonal(vars='All')
        # 
        # #nilonmin = np.argmax(dat.attr['lons']>loninf)-1
        # #nilonmax = np.argmax(dat.attr['lons']>=lonsup)+1
        # nilatmin = np.argmax(dat.attr['lats']>latinf)-1
        # #nilatmax = np.argmax(dat.attr['lats']>latsup)
        # if update_track==False:
        #     for altii in range(imaxalt,iminalt):
        #         for ilat in range(len(dats.var[variable][1])):
        #             for ilon in range(len(dats.var[variable][1][1])):
        #                 dats.var[variable][altii][ilat][ilon]=dats.var[variable][altii][ilat][ilon]-data_mean.var[variable][altii][nilatmin+ilat]
        # else:
        #     for ilat in range(len(dats.var[variable][1])):
        #         for ilon in range(len(dats.var[variable][1][1])): 
        #                 dats.var[variable][ialts][ilat][ilon]=dats.var[variable][ialts][ilat][ilon]-data_mean.var[variable][ialts][nilatmin+ilat]
        if update_track==False:
            min_O3=np.amin(dats.var[variable+'ano'][imaxalt:iminalt])
            maskmax=min_O3+min_O3*mask_percent
            maskmin=min_O3-min_O3*mask_percent
            argwhere=np.argwhere(dats.var[variable+'ano'][imaxalt:iminalt]==min_O3)    
            ialts=imaxalt+argwhere[0][0]
            ilats=argwhere[0][1]
            ilons=argwhere[0][2]
        ## Arrangement
        # else:
        #     ialts=np.argmax(dats.attr['zscale']<altitude)
        #     track_model[idx]['ialts']=ialts
        #     ilons=np.argmax(dats.attr['lons']>lon_mid)-1
        #     track_model[idx]['ilons']=ilons
        #     ilats=np.argmax(dats.attr['lats']>lat_bulle)-1
        #     track_model[idx]['ilats']=ilats
        #     min_O3=dats.var[variable+'ano'][ialts][ilats][ilons]
        #     track_model[idx]['min_O3']=min_O3
        print('min_O3 :',min_O3)
        
        # nilatmindats00 = np.argmax(85>latinf)-1
        # for ilat in range(len(dats00.var[variable][1])):
        #     for ilon in range(len(dats00.var[variable][1][1])):
        #         dats00.var[variable][ialts][ilat][ilon]=dats00.var[variable][ialts][ilat][ilon]-data_mean.var[variable][ialts][nilatmindats00+ilat]
        
        #dats.var[variable]=np.ma.masked_greater(dats.var[variable],maskmax)
        #dats.var[variable]=np.ma.masked_less(dats.var[variable],maskmin)
        
        ax=dats00.show(variable+'ano',ialts,cLines=False,projec=None,show=False,txt=str(variable)+' alt')
        #datsvisu.show(variable,ialts,cLines=False,projec=None,show=False)
        ax2=dats00.chartlonz(variable+'ano',lat=dats.attr['lats'][ilats],levs=(ialts-20,ialts+20),show=False,txt=str(variable)+' lat')
        #datsvisu.chartlonz(variable,lat=dats.attr['lats'][ilats],levs=(ialts-6,ialts+6),show=False)
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
                dat._get_var(variable) 
                # dat._get_var('T')
                # dat._get_var('U') 
                # dat._get_var('V')
                # dat._get_var('VO')
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
                    try:
                        if maskmin>maskmax:
                            mask1 = (dats.var[variable+'ano'][imaxalt:iminalt] < maskmin) & (dats.var[variable+'ano'][imaxalt:iminalt] > maskmax) 
                        elif maskmin<maskmax:
                            mask1 = (dats.var[variable+'ano'][imaxalt:iminalt] > maskmin) & (dats.var[variable+'ano'][imaxalt:iminalt] < maskmax) 
                        min_O3=np.amin(dats.var[variable+'ano'][imaxalt:iminalt][mask1])
                    except(TypeError,ValueError):
                        min_O3=np.amin(dats.var[variable+'ano'][imaxalt:iminalt])
                    print('min_O3 ;',min_O3)
                    argwhere=np.argwhere(dats.var[variable+'ano'][imaxalt:iminalt]==min_O3)    
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
                    idx+=1
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


def toggle_selector(event):
    if event.key in ['c','C']:
        print('close all')
        plt.close('all')

dirout='/data/hlestrel/fig_bulle_'+initials+'/'

# track_modelX={}
# nametrackmodelX='/home/hlestrel/stageM1-vortex/track_model_'+initials+'debutest'
# if os.path.isfile(nametrackmodelX):
#     with gzip.open(nametrackmodelX,'rb') as f:
#         track_modelX = pickle.load(f)
if visualization_track:
    nametrackmodel='/home/hlestrel/stageM1-vortex/track_model_'+initials+'_BulleB2_ref_tot'
    inputest=False
    import matplotlib.colors as colors
    import cartopy.crs as ccrs
    from cartopy import feature
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
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
    lon_visu=70
    lat_visu=45
    alt_visu=15
    if os.path.isfile(nametrackmodel) :
        with gzip.open(nametrackmodel,'rb') as f:
            track_model = pickle.load(f)
        print(track_model)
        idx=np.int(input('indice du dictionnaire   '))
        while idx!=999:
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
            O3mean = np.mean(dat.var['O3'],axis=2)
            Tmean = np.mean(dat.var['T'],axis=2)
            PVmean = np.mean(dat.var['PV'],axis=2)
            dat.var['PVano'] = dat.var['PV'] - PVmean[:,:,np.newaxis]
            dat.var['O3ano'] = dat.var['O3'] - O3mean[:,:,np.newaxis]
            dat.var['Tano'] = dat.var['T'] - Tmean[:,:,np.newaxis]
            #longrmin=-179
            #lonrngmx=179
            if (lonrngmin < 0):
                datp=dat.shift2west(180) 
                datsvisu= datp.extract(varss='All',lonRange=(lonrngmin,lonrngmx),latRange=(lat_bulle-lat_visu/2,lat_bulle+lat_visu/2))
                #datsvisu= datp.extract(varss='All',lonRange=(lonrngmin,lonrngmx),latRange=(10,88))
            else:
                datsvisu = dat.extract(varss='All',lonRange=(lonrngmin,lonrngmx),latRange=(lat_bulle-lat_visu/2,lat_bulle+lat_visu/2))
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
            ialts=np.argmax(dat.attr['zscale']<altitude)
            if inputest:
                ialts=np.argmax(dat.attr['zscale']<altitude)
            altmax=np.argmax(dat.attr['zscale']<altitude+alt_visu*0.3)
            if altitude-alt_visu*0.7>0:
                altmin=np.argmax(dat.attr['zscale']<altitude-alt_visu*0.7)
            else:
                altmin=np.argmax(dat.attr['zscale']<=1)
            ax=datsvisu.show('PVano',ialts,projec=None,show=False,txt='Anomaly of potential vorticity at altitude '+str(round(altitude,2)))#,clim=(-0.1*10**-5,1.5*10**-5))#,savfile=dirout+'PValt-'+str(track_model[idx]['date'])+str(idx)+'.png'
            # if lon_mid>180:
            #     ax.scatter(lon_mid+180,lat_bulle,marker='+',linewidth=4,color='yellow')
            #     #ax3.scatter(lon_mid+180,lat_bulle,marker='+',linewidth=4,color='yellow')
            #     #ax5.scatter(lon_mid+180,lat_bulle,marker='+',linewidth=4,color='yellow')
            # else:
            ax.scatter(lon_mid,lat_bulle,marker='+',linewidth=4,color='yellow')
            savfile=dirout+'PVanoalt-'+str(date00)+str(idx)+'.png'
            plt.savefig(savfile,bbox_inches='tight',dpi=300,format='png')
            ax2=datsvisu.chartlonz('PV',lat=lat_bulle,levs=(altmax,altmin),show=False,txt='Anomaly of potential vorticity')#,clim=(-0.1*10**-5,1.5*10**-5))#,savfile=dirout+'PVlat-'+str(track_model[idx]['date'])+str(idx)+'.png'
            if lon_mid>180:
                ax2.scatter(lon_mid-360,altitude,marker='o',linewidth=4,color='yellow')
            else:
                ax2.scatter(lon_mid,altitude,marker='o',linewidth=4,color='yellow')
            savfile2=dirout+'PVanolat-'+str(date00)+str(idx)+'.png'
            plt.savefig(savfile2,bbox_inches='tight',dpi=300,format='png')
            ax3=datsvisu.show('Tano',ialts,projec=None,show=False,txt='Tano alt '+str(altitude))#,savfile=dirout+'Tanoalt-'+str(track_model[idx]['date'])+str(idx)+'.png')#,clim=(-1,1))
            ax3.scatter(lon_mid,lat_bulle,marker='+',linewidth=4,color='yellow')
            savfile3=dirout+'Tanoalt-'+str(date00)+str(idx)+'.png'
            plt.savefig(savfile3,bbox_inches='tight',dpi=300,format='png')
            ax4=datsvisu.chartlonz('Tano',lat=lat_bulle,levs=(altmax,altmin),show=False,txt='Tano')#,savfile=dirout+'Tanolat-'+str(track_model[idx]['date'])+str(idx)+'.png')#,clim=(-1,1))
            if lon_mid>180:
                ax4.scatter(lon_mid-360,altitude,marker='o',linewidth=4,color='yellow')
            else:
                ax4.scatter(lon_mid,altitude,marker='o',linewidth=4,color='yellow')
            savfile4=dirout+'Tanolat-'+str(date00)+str(idx)+'.png'
            plt.savefig(savfile4,bbox_inches='tight',dpi=300,format='png')
            ax5=datsvisu.show('O3ano',ialts ,projec=None,show=False,txt=str('O3ano')+' alt '+str(altitude))#,savfile=dirout+'O3anoalt-'+str(track_model[idx]['date'])+str(idx)+'.png')#,clim=(-10*10**-7,-3*10**-7))
            ax5.scatter(lon_mid,lat_bulle,marker='+',linewidth=4,color='yellow')
            savfile5=dirout+'O3anoalt-'+str(date00)+str(idx)+'.png'
            plt.savefig(savfile5,bbox_inches='tight',dpi=300,format='png')
            ax6=datsvisu.chartlonz('O3ano',lat=lat_bulle,levs=(altmax,altmin),show=False,txt='O3ano')#,savfile=dirout+'O3anolat-'+str(track_model[idx]['date'])+str(idx)+'.png')#,clim=(-10*10**-7,-3*10**-7))
            if lon_mid>180:
                ax6.scatter(lon_mid-360,altitude,marker='o',linewidth=4,color='yellow')
            else:
                ax6.scatter(lon_mid,altitude,marker='o',linewidth=4,color='yellow')
            savfile6=dirout+'O3anolat-'+str(date00)+str(idx)+'.png'
            plt.savefig(savfile6,bbox_inches='tight',dpi=300,format='png')
            plt.connect('key_press_event', toggle_selector)
            plt.show()
            ## track_modelX[idx] = {'date':date00,'alt':altitude,'lon':lon_mid,'lat':lat_bulle,'min_O3':0,'ialts':0,'ilats':0,'ilons':0}
            ## with gzip.open(nametrackmodelX,'wb') as f:
            ##     pickle.dump(track_modelX,f)
            idx=np.int(input('indice du dictionnaire   '))
            
if print_traj:
    nametrackmodel1 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_A.pkl'
    nametrackmodel2 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_B1.pkl'
    nametrackmodel3 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_B2.pkl'
    savfile=diroutput
    variable='O3'
    fenetre_lon=6 # fenêtre DE REFERENCE de recherche de la bulle
    fenetre_lat=5
    fenetre_alt=2
    ###### UPDATE TRACK
    colori=170
    track_model1 = pickle.load(open(nametrackmodel1,'rb'))
    track_model2 = pickle.load(open(nametrackmodel2,'rb'))
    track_model3 = pickle.load(open(nametrackmodel3,'rb'))
    idx1=len(track_model1['dates'])
    idx2=len(track_model2['dates'])
    idx3=len(track_model3['dates'])
    date00=track_model2['dates'][200] #choix de l'indice au niveau de la séparation en trois du panache
    altitude=track_model2['alts'][200]
    lon_mid=track_model2['lons'][200]
    lat_bulle=track_model2['lats'][200]
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
    if os.path.isfile(savfile+'dats00traj_'+initials):
        with gzip.open(savfile+'dats00traj_'+initials,'rb') as f:
            dats00=pickle.load(f)
        ialts=np.argmax(dats00.attr['zscale']<altitude)
    else:
        dat = ECMWF('FULL-EA',date00,exp='VOZ') 
        dat._get_var(variable)
        dat._mkpscale() 
        dat._mkzscale()
        ialts=np.argmax(dat.attr['zscale']<altitude)
        O3mean = np.mean(dat.var[variable],axis=2)
        dat.var[variable+'ano'] = dat.var[variable] - O3mean[:,:,np.newaxis]
        #if (loninf < 0):
        datp=dat.shift2west(180) 
        dats = datp.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
        if lat_bulle>0:
            dats00 = datp.extract(varss='All',lonRange=(-180,179),latRange=(10,85))
        elif lat_bulle<0:
            dats00 = datp.extract(varss='All',lonRange=(-180,179),latRange=(-80,-3))
        with gzip.open(savfile+'dats00traj_'+initials,'wb') as f:
            pickle.dump(dats00,f)
    
    #datsvisu.chartlonz(variable,lat=dats.attr['lats'][ilats],levs=(ialts-6,ialts+6),show=False)
    jet1= plt.get_cmap('Reds') # jet
    colors1 = iter(jet1(np.linspace(0.1,1,colori)))
    colors12 = iter(jet1(np.linspace(0.1,1,colori))) # attention à bien adapter au nombre de points que l'on souhaite tracer
    jet2= plt.get_cmap('Blues') 
    colors2 = iter(jet2(np.linspace(0.1,1,colori))) 
    colors22 = iter(jet2(np.linspace(0.1,1,colori))) 
    jet3= plt.get_cmap('binary')
    colors3 = iter(jet3(np.linspace(0.1,1,colori))) 
    colors32 = iter(jet3(np.linspace(0.1,1,colori))) 
    lons1=[]
    lats1=[]
    alts1=[]
    for i in range(idx1):
        alts1.append(track_model1['PT'][i])
        lats1.append(track_model1['lats'][i])
        if track_model1['lons'][i]>180:
            lons1.append(track_model1['lons'][i]-360)
        elif track_model1['lons'][i]<-180:
            lons1.append(track_model1['lons'][i]+360)
        else:
            lons1.append(track_model1['lons'][i])
    lons2=[]
    lats2=[]
    alts2=[]
    for i in range(idx2):
        alts2.append(track_model2['PT'][i])
        lats2.append(track_model2['lats'][i])
        if track_model2['lons'][i]>180:
            lons2.append(track_model2['lons'][i]-360)
        elif track_model2['lons'][i]<-180:
            lons2.append(track_model2['lons'][i]+360)
        else:
            lons2.append(track_model2['lons'][i])
    lons3=[]
    lats3=[]
    alts3=[]
    for i in range(idx3):
        alts3.append(track_model3['PT'][i])
        lats3.append(track_model3['lats'][i])
        if track_model3['lons'][i]>180:
            lons3.append(track_model3['lons'][i]-360)
        elif track_model3['lons'][i]<-180:
            lons3.append(track_model3['lons'][i]+360)
        else:
            lons3.append(track_model3['lons'][i])
    itera=2
    lnwd=3
    
##    ax2=dats00.chartlonz(variable+'ano',lat=dats.attr['lats'][ilats],levs=(ialts-20,ialts+20),show=False,txt=str(variable)+' lat')
    lat=lat_bulle
    figsize=(11,4)
    #fig = plt.figure(figsize=figsize,constrained_layout=True)
    levs=(ialts-20,ialts+20)
    #txt='Vortex trajectories zonal section'
    log=False
    clim=(None,None)
    cmap=mymap
    scale=10**6
    try:
        pos=np.where(dats00.attr['lats']>=lat)[0][0]
        print('pos',pos)
    except:
        print('lat out of range')
    fig2 = plt.figure(figsize=figsize)
    fs = 15
    ax2 = fig2.add_subplot(111)
    if levs[0]==None: l1=29
    else: l1 = levs[0]
    if levs[1]==None: l2=115
    else: l2 = levs[1]
    lons = np.arange(dats00.attr['lons'][0]-0.5*dats00.attr['dlo'],dats00.attr['lons'][-1]+dats00.attr['dlo'],dats00.attr['dlo'])
    # try:
    #     zz1 = 0.5*(dats00.var['Z'][l1-1:l2+1,pos, :] + dats00.var['Z'][l1:l2+2,pos,:])/1000
    #     zz = np.empty((zz1.shape[0],zz1.shape[1]+1))
    #     zz[:,1:-1] = 0.5*(zz1[:,1:]+zz1[:,:-1])
    #     zz[:,0] = zz1[:,0]
    #     zz[:,-1] = zz1[:,-1]
    #     print(zz.shape,len(lons))
    #     iax2=ax2.pcolormesh(lons,zz,scale*dats00.var[variable+'ano'][l1:l2+1,pos,:],
    #                     vmin=clim[0],vmax=clim[1],cmap=cmap)
    #     plt.ylabel('altitude (km)',fontsize=fs)
    #     print('USE Z')
    # except(KeyError):
    #     iax2=ax2.pcolormesh(lons,dats00.attr['zscale_i'][l1:l2+2],dats00.var[variable+'ano'][l1:l2+1,pos, :],
    #                 vmin=clim[0],vmax=clim[1],cmap=cmap)

    #       plt.ylabel('baro altitude (km)',fontsize=fs)
    plt.ylabel('potential temperature (K)',fontsize=fs)
    ax2.tick_params(labelsize=16)
    plt.xlabel('longitude',fontsize=fs)

    #plt.title(txt+" latitude "+str(lat)+'°',fontsize=fs)
    #cbar = fig2.colorbar(iax2)
    #cbar.ax.tick_params(labelsize=fs)
    
    for i in range(itera,idx1,itera):
        #if  54<=i<=248:
        if  track_model2['dates'][i].month==9 and 3<track_model2['dates'][i].day<16:
            #c=next(colors12)
            if (lons1[i-itera]-lons1[i]<-100): 
                ax2.plot([180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c,marker='x')
            elif (lons1[i-itera]>150 and lons1[i])<-150:
                ax2.plot([-180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c,marker='x')
            #elif i==200:
                #ax2.plot([lons1[200-itera],lons1[200]],[alts1[200-itera],alts1[200]],linewidth=lnwd,color=c,marker='x',label='no CALIOP')
            else:
                ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c,marker='x')
        elif i==int(idx1/itera)*itera:
            c=next(colors12)
            if (lons1[i-itera]-lons1[i]<-100): 
                ax2.plot([180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
            elif (lons1[i-itera]>150 and lons1[i])<-150:
                ax2.plot([-180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
            else:
                ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
        else:
            c=next(colors12)
            if (lons1[i-itera]-lons1[i]<-100): 
                ax2.plot([180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
            elif (lons1[i-itera]>150 and lons1[i])<-150:
                ax2.plot([-180,lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
            else:
                ax2.plot([lons1[i-itera],lons1[i]],[alts1[i-itera],alts1[i]],linewidth=lnwd,color=c)
    ax2.plot([lons1[len(lons1)-1-itera],lons1[len(lons1)-1]],[alts1[len(lons1)-1-itera],alts1[len(lons1)-1]],linewidth=lnwd,color=c,label='Vortex A')
    for i in range(itera,idx2,itera):
        if  track_model2['dates'][i].month==9 and 3<track_model2['dates'][i].day<16:
            if (lons2[i-itera]-lons2[i]<-100): 
                ax2.plot([180,lons2[i]],[alts2[i-itera],alts2[i]],linewidth=lnwd,color=c,marker='x')
            elif (lons2[i-itera]>150 and lons2[i])<-150:
                ax2.plot([-180,lons2[i]],[alts2[i-itera],alts2[i]],linewidth=lnwd,color=c,marker='x')
            else:
                ax2.plot([lons2[i-itera],lons2[i]],[alts2[i-itera],alts2[i]],linewidth=lnwd,color=c,marker='x')
        elif i==int(idx2/itera)*itera:
            c=next(colors22)
            if (lons2[i-itera]-lons2[i-1]<-100): 
                ax2.plot([180,lons2[i-1]],[alts2[i-itera],alts2[i-1]],linewidth=lnwd,color=c)
            elif (lons2[i-itera]>150 and lons2[i-1])<-150:
                ax2.plot([-180,lons2[i-1]],[alts2[i-itera],alts2[i-1]],linewidth=lnwd,color=c)
            else:
                ax2.plot([lons2[i-itera],lons2[i-1]],[alts2[i-itera],alts2[i-1]],linewidth=lnwd,color=next(colors22))
        else:
            c=next(colors22)
            if (lons2[i-itera]-lons2[i]<-100): 
                ax2.plot([180,lons2[i]],[alts2[i-itera],alts2[i]],linewidth=lnwd,color=c)
            elif (lons2[i-itera]>150 and lons2[i])<-150:
                ax2.plot([-180,lons2[i]],[alts2[i-itera],alts2[i]],linewidth=lnwd,color=c)
            else:
                ax2.plot([lons2[i-itera],lons2[i]],[alts2[i-itera],alts2[i]],linewidth=lnwd,color=c)
    ax2.plot([lons2[len(lons2)-1-itera],lons2[len(lons2)-2]],[alts2[len(lons2)-itera-1],alts2[len(lons2)-2]],linewidth=lnwd,color=c,label='Vortex B1')

    for i in range(itera,idx3,itera):
        if  track_model3['dates'][i].month==9 and 3<track_model3['dates'][i].day<16:
            if (lons3[i-itera]-lons3[i]<-100): 
                ax2.plot([180,lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color=c,marker='x')
            elif (lons3[i-itera]>150 and lons3[i])<-150:
                ax2.plot([-180,lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color=c,marker='x')
            elif i==200:
                    ax2.plot([lons3[200-itera],lons3[200]],[alts3[200-itera],alts3[200]],linewidth=lnwd,color=c,marker='x',label='no CALIOP')
            else:
                ax2.plot([lons3[i-itera],lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color=c,marker='x')
        elif i==int(idx3/itera)*itera:
            c=next(colors32)
            if (lons3[i-itera]-lons3[i]<-100): 
                ax2.plot([180,lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color=c)
            elif (lons3[i-itera]>150 and lons3[i])<-150:
                ax2.plot([-180,lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color=c)
            else:
                ax2.plot([lons3[i-itera],lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color=c)
        else:
            c=next(colors32)
            if (lons3[i-itera]-lons3[i]<-100): 
                ax2.plot([180,lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color=c)
            elif (lons3[i-itera]>150 and lons3[i])<-150:
                ax2.plot([-180,lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color=c)
            else:
                ax2.plot([lons3[i-itera],lons3[i]],[alts3[i-itera],alts3[i]],linewidth=lnwd,color=c)
    ax2.plot([lons3[len(lons3)-1-itera],lons3[len(lons3)-1]],[alts3[len(lons3)-itera-1],alts3[len(lons3)-1]],linewidth=lnwd,color=c,label='Vortex B2')

    plt.xticks([-180,-120,-60,0,60,120,180],['180°W','120°W','60°W','0°','60°E','120°E','180°E'])
    ax2.set_xlim(-180,180)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,fontsize=fs)
    plt.savefig(savfile+str(variable)+'_lat_PT_'+initials+'.png',bbox_inches='tight',dpi=300,format='png')
##
    
    ##ax=dats00.show(variable+'ano',ialts,cLines=False,projec=None,show=False,txt=str(variable)+' alt')
    lev=ialts
    cardinal_level=True
    #txt='Vortex trajectories horizontal section at '+str(round(altitude,0))+' km'
    log=False
    clim=(None,None)
    cmap=mymap
    cLines=None
    scale=10**6
    aspect=1
    sat_H=35785831
    xylim=False
    if len(dats00.var[variable+'ano'].shape) == 3:
        if (cardinal_level==False) | (lev > dats00.nlev-1):
            clev = np.abs(np.array(dats00.attr['levs'])-lev).argmin()
        else:
            clev = lev
        buf = dats00.var[variable+'ano'][clev,:,:]
    else:
        buf = dats00.var[variable+'ano']
    fs=15
    cm_lon =0
    ctrl_lon=(dats00.attr['lons'][0]+dats00.attr['lons'][-1])/2
    ctrl_lat=(dats00.attr['lats'][0]+dats00.attr['lats'][-1])/2
    if dats00.attr['lons'][-1] > 180: cm_lon=180
    projplate = ccrs.PlateCarree(central_longitude=cm_lon)
    proj=projplate
    fig1 = plt.figure(figsize=(11,4))
    fig1.subplots_adjust(hspace=0,wspace=0.5,top=0.925,left=0.)
    ax = plt.axes(projection = proj)
    #iax = ax.imshow(scale*buf, transform=projplate,cmap=cmap,clim=clim,extent=[dats00.attr['lons'][0]-cm_lon, dats00.attr['lons'][-1]-cm_lon,dats00.attr['lats'][0], dats00.attr['lats'][-1]],interpolation='nearest',origin='lower',aspect=aspect)
    if xylim:
        x1a,y1a = proj.transform_point(dats00.attr['lons'][0]-2*cm_lon,dats00.attr['lats'][0],ccrs.Geodetic())
        x1b,y1b = proj.transform_point(dats00.attr['lons'][0]-2*cm_lon,dats00.attr['lats'][-1],ccrs.Geodetic())
        x1c,y1c = proj.transform_point(ctrl_lon-2*cm_lon,dats00.attr['lats'][0],ccrs.Geodetic())
        x1 = min(x1a,x1b,x1c)
        y1 = min(y1a,y1b,y1c)
        x2a,y2a = proj.transform_point(dats00.attr['lons'][-1]-2*cm_lon,dats00.attr['lats'][-1],ccrs.Geodetic())
        x2b,y2b = proj.transform_point(dats00.attr['lons'][-1]-2*cm_lon,dats00.attr['lats'][0],ccrs.Geodetic())
        x2c,y2c = proj.transform_point(ctrl_lon-2*cm_lon,dats00.attr['lats'][-1],ccrs.Geodetic())
        x2 = max(x2a,x2b,x2c)
        y2 = max(y2a,y2b,y2c)
        ax.set_xlim(x1,x2)
        ax.set_ylim(y1,y2)
    xlocs = None
    if cm_lon == 180:
            interx = 30
            minx = dats00.attr['lons'][0] + interx - dats00.attr['lons'][0]%interx
            xlocs = list(np.arange(minx,181,interx))+list(np.arange(interx-180,self.attr['lons'][-1]-360,interx))
    gl = ax.gridlines(draw_labels=True, xlocs=xlocs,
                    linewidth=2, color='gray', alpha=0.5, linestyle='--')
    if cLines is not None:
                ax.contour(dats00.var[variable+'ano'][lev,:,:],transform=projplate,extent=(dats00.attr['lons'][0]-cm_lon,
                            dats00.attr['lons'][-1]-cm_lon,dats00.attr['lats'][0],dats00.attr['lats'][-1]),levels=cLines,origin='lower')
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
    #plt.title(txt,fontsize=fs)
    # axpos = ax.get_position()
    # pos_x = axpos.x0 + axpos.x0 + axpos.width + 0.01
    # pos_cax = fig1.add_axes([pos_x,axpos.y0,0.04,axpos.height])
    # cbar=fig1.colorbar(iax,cax=pos_cax)
    # cbar.ax.tick_params(labelsize=fs)
    ##
    
    for i in range(itera,idx1,itera):
        if  track_model1['dates'][i].month==9 and 3<track_model1['dates'][i].day<16:
        #if  54<=i<=248:
            #c=next(colors1)
            if (lons1[i-itera]-lons1[i]<-100): 
                ax.plot([180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x')
            elif (lons1[i-itera]>150 and lons1[i])<-150:
                ax.plot([-180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x')
            else:
                ax.plot([lons1[i-itera],lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c,marker='x')
        else:
            c=next(colors1)
            if (lons1[i-itera]-lons1[i]<-100): 
                ax.plot([180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
            elif (lons1[i-itera]>150 and lons1[i])<-150:
                ax.plot([-180,lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
            else:
                ax.plot([lons1[i-itera],lons1[i]],[lats1[i-itera],lats1[i]],linewidth=lnwd,color=c)
    for i in range(itera,idx2,itera):
        if  track_model2['dates'][i].month==9 and 3<track_model2['dates'][i].day<16:
            if (lons2[i-itera]-lons2[i]<-100): 
                ax.plot([180,lons2[i]],[lats2[i-itera],lats2[i]],linewidth=lnwd,color=c,marker='x')
            elif (lons2[i-itera]>150 and lons2[i])<-150:
                ax.plot([-180,lons2[i]],[lats2[i-itera],lats2[i]],linewidth=lnwd,color=c,marker='x')
            else:
                ax.plot([lons2[i-itera],lons2[i]],[lats2[i-itera],lats2[i]],linewidth=lnwd,color=c,marker='x')
        else:
            c=next(colors2)
            if (lons2[i-itera]-lons2[i]<-100): 
                ax.plot([180,lons2[i]],[lats2[i-itera],lats2[i]],linewidth=lnwd,color=c)
            elif (lons2[i-itera]>150 and lons2[i])<-150:
                ax.plot([-180,lons2[i]],[lats2[i-itera],lats2[i]],linewidth=lnwd,color=c)
            else:
                ax.plot([lons2[i-itera],lons2[i]],[lats2[i-itera],lats2[i]],linewidth=lnwd,color=c)
    for i in range(itera,idx3,itera):
        if  track_model3['dates'][i].month==9 and 3<track_model3['dates'][i].day<16:
            if (lons3[i-itera]-lons3[i]<-100): 
                ax.plot([180,lons3[i]],[lats3[i-itera],lats3[i]],linewidth=lnwd,color=c,marker='x')
            elif (lons3[i-itera]>150 and lons3[i])<-150:
                ax.plot([-180,lons3[i]],[lats3[i-itera],lats3[i]],linewidth=lnwd,color=c,marker='x')
            else:
                ax.plot([lons3[i-itera],lons3[i]],[lats3[i-itera],lats3[i]],linewidth=lnwd,color=c,marker='x')
        else:
            c=next(colors3)
            if (lons3[i-itera]-lons3[i]<-100): 
                ax.plot([180,lons3[i]],[lats3[i-itera],lats3[i]],linewidth=lnwd,color=c)
            elif (lons3[i-itera]>150 and lons3[i])<-150:
                ax.plot([-180,lons3[i]],[lats3[i-itera],lats3[i]],linewidth=lnwd,color=c)
            else:
                ax.plot([lons3[i-itera],lons3[i]],[lats3[i-itera],lats3[i]],linewidth=lnwd,color=c)
    plt.savefig(savfile+str(variable)+'_horiz_PT_'+initials+'.png',dpi=300,bbox_inches='tight',format='png')
    plt.show()



# nametrackmodel1 = '/home/hlestrel/stageM1-vortex/track_model_AUS2009_mean_Hugo'
# with gzip.open(nametrackmodel1,'rb') as f:
#     track_model1 = pickle.load(f)
#     
# track_modelY={}
# track_modelY['dates']=[]
# track_modelY['lons']=[]
# track_modelY['lats']=[]
# track_modelY['alts']=[]
# for i in range(len(track_model1)):
#     track_modelY['dates'].append(track_model1[i]['date'])
#     track_modelY['lons'].append(track_model1[i]['lon'])
#     track_modelY['lats'].append(track_model1[i]['lat'])
#     track_modelY['alts'].append(track_model1[i]['alt'])
# pickle.dump(track_modelY,open('/home/hlestrel/stageM1-vortex/Vortex-track_AUS2009.pkl','wb'))
# 
# initials='CAN2017'
# nametrackmodel1 = '/home/hlestrel/stageM1-vortex/Vortex-track_'+initials+'_A.pkl'
# track_model1 = pickle.load(open(nametrackmodel1,'rb'))
# date00=track_model1['dates'][-1]
# dat = ECMWF('FULL-EA',date00,exp='VOZ') 
# dat._get_var('O3')
# dat._get_var('T')
# dat._get_var('U') 
# dat._get_var('V')
# dat._get_var('VO')
# dat._get_var('P')
# dat._get_var('PT')
# dat._mkpscale() 
# dat._mkzscale()
# dat._mkp() 
# dat._mkthet() 
# dat._mkpv()
# PTmean = np.mean(dat.var['PT'],axis=2)
# Pmean2=np.mean(dat.var['P'],axis=(1,2))
# dat.var['PTano'] = dat.var['PT'] - PTmean[:,:,np.newaxis]
# dats = dat.extract(varss='All',lonRange=(0,359),latRange=(-85,85)) 
# dats.attr['pzscale']=np.log(Pmean/100000)*(-7.4)  
# 
# track_model1['PT']=[]
# for i in range(len(track_model1['lats'])):
#     track_model1['PT'].append(PTmean[np.where(dats.attr['zscale']<=track_model1['alts'][i])[0][0]][np.where(dats.attr['lats']>=track_model1['lats'][i])[0][0]])
#     
# pickle.dump(track_model1,open(nametrackmodel1,'wb'))


# date00=datetime(2017,9,2,3)
# #date00=track_model1['dates'][-1]
# dat = ECMWF('FULL-EA',date00,exp='VOZ') 
# dat._get_var('O3')
# dat._get_var('T')
# dat._get_var('U') 
# dat._get_var('V')
# dat._get_var('VO')
# dat._get_var('P')
# dat._get_var('PT')
# dat._mkpscale() 
# dat._mkzscale()
# dat._mkp() 
# dat._mkthet() 
# dat._mkpv()
# PVmean = np.mean(dat.var['PV'],axis=2)
# Tmean = np.mean(dat.var['T'],axis=2)
# Pmean=np.mean(dat.var['P'],axis=(1,2))
# dat.var['PVano'] = dat.var['PV'] - PVmean[:,:,np.newaxis]
# dat.var['Tano'] = dat.var['T'] - Tmean[:,:,np.newaxis]
# #datp=dat.shift2west(210) 
# #dats = datp.extract(varss='All',lonRange=(-150,209),latRange=(-25,89))
# dats = dat.extract(varss='All',lonRange=(0,359),latRange=(-25,89))
# #dats.attr['pzscale']=np.log(Pmean/100000)*(-7.4)   
# ial=np.where(dats.attr['zscale']<=19)[0][0]
# #mymap=plt.cm.get_cmap('viridis')#, 10000)#'RdGy'
# dats.show('PVano',ial,projec='nearside',sat_H=9578583,scale=10**5,clim=(-0.8,0),cmap=plt.cm.get_cmap('bone'),figsize=(22,8),savfile='image_page_rapport.png')
# #dats.chartlonz('Tano',lat=50,levs=(ial-5,ial+10))