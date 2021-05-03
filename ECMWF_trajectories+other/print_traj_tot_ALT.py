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

initials='CAN2017'
diroutput='/data/hlestrel/fig_bulle_'+initials+'/'

# ERA5 track
track_model1 = pickle.load(open(diroutput+'print_traj_tot/Vortex-track_'+initials+'_A.pkl','rb'))
track_model2 = pickle.load(open(diroutput+'print_traj_tot/Vortex-track_'+initials+'_B1.pkl','rb'))
track_model3 = pickle.load(open(diroutput+'print_traj_tot/Vortex-track_'+initials+'_B2.pkl','rb'))
# CALIOP track
with gzip.open(diroutput+'track_CALIOP/Ntrack-L1-'+initials,'rb') as f:
    CALIOPtrack=pickle.load(f)
events=[65,79,93,107,149,178,187,201,215,216,229,230,231,244,259,274,288,317,331,358,411,425,439,453,593,607,664,678,692]
Ctracklon=[]
Ctracklat=[]
Ctrackalt=[]
#for i in range(len(CALIOPtrack)):
for j in events:
    for i in range(len(CALIOPtrack)):
        if CALIOPtrack[i]['cald_index']==j:
            break
    Ctrackalt.append(CALIOPtrack[i]['mid'])
    Ctracklat.append(CALIOPtrack[i]['lat_mid'])
    if CALIOPtrack[i]['lon_mid']>180:
        Ctracklon.append(CALIOPtrack[i]['lon_mid']-360)
    elif CALIOPtrack[i]['lon_mid']<-180:
        Ctracklon.append(CALIOPtrack[i]['lon_mid']+360)
    else:
        Ctracklon.append(CALIOPtrack[i]['lon_mid'])
savfile=diroutput+'print_traj_tot/'
variable='O3'
fenetre_lon=6 # fenêtre DE REFERENCE de recherche de la bulle
fenetre_lat=5
fenetre_alt=2
colori=170
idx1=len(track_model1['dates'])
idx2=len(track_model2['dates'])
idx3=len(track_model3['dates'])
# date00=track_model2['dates'][200] #choix de l'indice au niveau de la séparation en trois du panache
# altitude=track_model2['alts'][200]
# lon_mid=track_model2['lons'][200]
# lat_bulle=track_model2['lats'][200]
date00=track_model1['dates'][67] #choix de l'indice au niveau de la séparation en trois du panache
altitude=track_model1['alts'][67]
lon_mid=track_model1['lons'][67]
lat_bulle=track_model1['lats'][67]
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
    dat._get_var('O3')
    dat._get_var('Z')
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
    Pmean=np.mean(dat.var['P'],axis=(1,2))
    ialts=np.argmax(dat.attr['zscale']<altitude)
    O3mean = np.mean(dat.var['O3'],axis=2)
    dat.var['O3ano'] = dat.var['O3'] - O3mean[:,:,np.newaxis]
    PTmean = np.mean(dat.var['PT'],axis=2)
    dat.var['PTmean']=dat.var['PT']-PTmean[:,:,np.newaxis]
    #if (loninf < 0):
    datp=dat.shift2west(180) 
    dats = datp.extract(varss='All',lonRange=(loninf,lonsup),latRange=(latinf,latsup)) 
    if lat_bulle>0:
        dats00 = datp.extract(varss='All',lonRange=(-180,179),latRange=(10,85))
    elif lat_bulle<0:
        dats00 = datp.extract(varss='All',lonRange=(-180,179),latRange=(-80,-3))
    dats00.var['PTmean']=np.mean(dats.var['PTmean']+dats.var['PT'],axis=2)
    dats00.attr['pzscale']=np.log(Pmean/100000)*(-7.4)
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
    alts1.append(dats00.attr['pzscale'][np.where(dats00.attr['zscale']<=track_model1['alts'][i])[0][0]])
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
    alts2.append(dats00.attr['pzscale'][np.where(dats00.attr['zscale']<=track_model2['alts'][i])[0][0]])
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
    alts3.append(dats00.attr['pzscale'][np.where(dats00.attr['zscale']<=track_model3['alts'][i])[0][0]])
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
plt.ylabel('Altitude (km)',fontsize=fs)
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

ax2.scatter(Ctracklon,Ctrackalt,linewidth=2,color='g',marker='x')

plt.xticks([-180,-120,-60,0,60,120,180],['180°W','120°W','60°W','0°','60°E','120°E','180°E'])
ax2.set_xlim(-180,180)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,fontsize=fs)
plt.savefig(savfile+str(variable)+'_lon_PT_'+initials+'.png',bbox_inches='tight',dpi=300,format='png')
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

#PLOT CALIOP TRACK
ax.scatter(Ctracklon,Ctracklat,linewidth=2,color='g',marker='x')

plt.savefig(savfile+str(variable)+'_horiz_PT_'+initials+'.png',dpi=300,bbox_inches='tight',format='png')
plt.show()
