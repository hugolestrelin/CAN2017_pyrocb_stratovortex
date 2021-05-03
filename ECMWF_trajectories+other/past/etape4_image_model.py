#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Shows the model from the data of the chosen sequences of CALIOP sections in function of the track dictionary (trackname), from the directory defined in diroutput.
This opration is realised in a separated serveur, the tracks have therefore to be previously copied from a serveur to another. 

input : track dictionary of chosen sequences (diroutput ; trackname)
output : plots of the model in diroutput

@author: Bernard Legras, modified by Hugo Lestrelin

"""
import sys,os  
import socket      
if 'ciclad' in socket.gethostname():
    import sys,os
    sys.path.append('/home/hlestrel/pylib')                
from datetime import datetime,timedelta
from ECMWF_N_test_proj import ECMWF 
import numpy as np 
import matplotlib.pyplot as plt 
import pickle,gzip
from matplotlib import gridspec
import matplotlib.gridspec as gridspec
from progress_bar import printProgressBar

#### AJOUTER COMPOSITE ET TRACE DE LA SELCTION (e.g point en mid)
figsave=False 
show=False
mean=True
figsavemean=True 
showmean=False
composite=False
initials='AUS2009'
if initials=='AUS2009':
    datefirst=datetime(2009,1,28)
diroutput='/data/hlestrel/fig_bulle_'+initials+'/'
trackname='Ntrack-L1-'+initials
with gzip.open(diroutput+trackname,'rb') as f:
    infos = pickle.load(f)
# size of the lon/lat window
fenetre_lon=55 
fenetre_lat=40

# Does the update is needed or not :
update_index=[]
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
        date=dd
    
    
    lon_mid=infos[i]['lon_mid']
    #if lon_mid<(180-fenetre_lon/2):
    loninf = lon_mid-fenetre_lon/2
    lonsup = lon_mid+fenetre_lon/2
    #else:
        #loninf = (lon_mid-(360-fenetre_lon/2))-fenetre_lon/2
        #lonsup = (lon_mid-(360-fenetre_lon/2))+fenetre_lon/2
        
    lat_bulle =infos[i]['lat_mid']
    if (lat_bulle-fenetre_lat/2)>-85 or (lat_bulle+fenetre_lat/2)<85:
        latinf = lat_bulle-fenetre_lat/2
        latsup = lat_bulle+fenetre_lat/2
    elif (lat_bulle-fenetre_lat/2)<-85:
        latinf = -85
        latsup = lat_bulle+fenetre_lat/2
    elif (lat_bulle+fenetre_lat/2)>85:
        latinf = lat_bulle-fenetre_lat/2
        latsup = 85
    
    altitude = infos[i]['mid']
    
    dat = ECMWF('FULL-EA',date,exp='VOZ') 
    dat._get_var('U') 
    dat._get_var('V')
    dat._get_var('T')
    dat._get_var('VO') 
    dat._get_var('P')
    dat._get_var('PT') 
    try:
        dat._mkpscale() 
        dat._mkzscale()
        dat._mkp() 
        dat._mkthet() 
        dat._mkpv()
        ialt=0
        while dat.attr['zscale'][ialt]>altitude:
            ialt+=1
        if (loninf < 0):
            datp=dat.shift2west(150) 
            dats = datp.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
        else:
            dats = dat.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
            
        dats.var['PV']=dats.var['PV']*10**6
        
        if figsave:
            dats.show('PV',ialt,cLines=False,projec='nearside',show=False,savfile=diroutput+'plotVPalt-'+str(infos[i]['cald_index'])+'.png',txt='Potential vorticity alt = '+str(altitude)+'km '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
            plt.close('all')
            dats.chartlonz('PV',lat=lat_bulle,levs=(ialt-10,ialt+10),show=False,savfile=diroutput+'plotVPlat-'+str(infos[i]['cald_index'])+'.png',txt='Potential vorticity '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
            plt.close('all')
            dats.chartlonz('T',lat=lat_bulle,levs=(ialt-10,ialt+10),show=False,savfile=diroutput+'plotTlat-'+str(infos[i]['cald_index'])+'.png',txt='Temperature '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
            plt.close('all')
            print('The plot for the index '+str(infos[i]['cald_index'])+' has been saved')
        if show:
            dats.show('PV',ialt,cLines=False,projec='nearside',txt='Potential vorticity alt = %3.2f km '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour) % (altitude))
            dats.chartlonz('PV',lat=lat_bulle,levs=(ialt-10,ialt+10),txt='Potential vorticity '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
            dats.chartlonz('T',lat=lat_bulle,levs=(ialt-10,ialt+10),txt='Temperature '+str(date.year)+'-'+str(date.month)+'-'+str(date.day)+'-'+str(date.hour))
    
        if mean:
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
    except('U not found or read error'):
        print('Problème avec les données de l"index :',infos[i]['cald_index'])
        continue
        
#### AJOUTER LA FONCTION COMPOSITE :
    if composite:
        a=1
        