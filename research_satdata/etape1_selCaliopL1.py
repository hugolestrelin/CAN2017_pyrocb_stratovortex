#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Exploring the calipso database for orbits defined by the user in the box. First step of the finding_vortex program, needed in order to realise the plot in step 2.

input : nothing
output : dictonary object Cald named in filecal, in the directory defined by dirouput.

@author: Bernard Legras, modified by Hugo Lestrelin

"""
import numpy as np
from datetime import date,datetime,timedelta
import pickle,gzip
import matplotlib.pyplot as plt
import os
from pyhdf.SD import SD
import glob
from astropy.time import Time, TimeDelta

initials='++CAN2017' # project initials, please always keep the year as the 4 last char, like : AUS2017.

# Main directories for aerosol profiles and L1 data on the ICARE database
if int(initials[-4:])>2019:
    dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v3.40' # Watch out, has to be constantly updated.
    dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1.v3.40'
else:
    dirAProf = '/DATA/LIENS/CALIOP/05kmAPro.v4.10'
    dirL1 = '/DATA/LIENS/CALIOP/CAL_LID_L1.v4.10'
   
diroutput = '/scratch/hlestrelin/fig_'+initials+'/'

# Here you can add some orbits boxes to be used again if needed
## Australia 2019
if initials=='AUS2019':
    box = [0,359.99,-85,0]
    date0 = datetime(2019,1,8)
    nbday = 1
## Canada 2017
elif initials=='CAN2017':
    box = [0,359.99,0,87]
    date0 = datetime(2017,10,16)
    nbday = 85
## Canada 2017 16/09-15/10
elif initials=='+CAN2017':
    box = [0,359.99,0,87]
    date0 = datetime(2017,9,16)
    nbday = 29
    #events=[]
## Canada 2017 30/08 - 05/09
elif initials=='++CAN2017':
    box = [0,359.99,0,87]
    date0 = datetime(2017,8,30)
    nbday = 7
    #events=[]
## Asutralie 2015
elif initials=='AUS2015':
    box = [0,359.99,-85,0]
    date0 = datetime(2015,1,28)
    nbday = 40
## Asutralie 2014
elif initials=='AUS2014':
    box = [0,359.99,-85,0]
    date0 = datetime(2014,1,15)
    nbday = 40
## Australia 2009
elif initials=='AUS2009':
    box = [0,359.99,-85,0]
    date0 = datetime(2009,1,28)
    nbday = 1
## Alaska 2014
elif initials=='ALA2014':
    box = [0,359.99,0,85]
    date0 = datetime(2014,7,2)
    nbday = 40
## Canada 2011
elif initials=='CAN2011':
    box = [0,359.99,0,85]
    date0 = datetime(2011,5,15)
    nbday = 40
## AUS2020
if initials=='AUS2020':
	box= [100,300,-55,-15]
	date0=datetime(2020,2,2)
	nbday=10	

# Does the update is need or not :
# Attention lors de la réécriture du dictionnaire de ne pas faire à nouveau les mêmes jours
# Ajouter une option qui l'empêche
filecal='selCaliop_'+initials+'.pkl'
namecal = os.path.join( diroutput, filecal )
if os.path.isfile(namecal):
    update = True
    print('Update = True')
    with gzip.open(namecal,'rb') as f:
        box,idate0,nbday,Cald = pickle.load(f)
    idx = len(Cald)
    print ('read delection with ',idx,' records')
    print ('indexing from ',idx)
    nbday=nbday-(date(Cald[idx]['utc'].year,Cald[idx]['utc'].month,Cald[idx]['utc'].day)-date(date0.year,date0.month,date0.day)).days
    date0=Cald[idx]['utc']
    date00=Cald[1]['utc']
else:
    # Creation of the Cald variable, That's where the informations access will be stored
    update = False
    print('Update = False, new dictionary')
    Cald = {}
    date00 = date0
    nbday0 = 0
    idx = 0

nbday0 = 0
print('date0,nbday,box')
print(date0,nbday,box)

# Browse dates
print(nbday,' days to be processed')
for iday in range(nbday):
    date = date0 + timedelta(days=iday)
    #if initials=='CAN2017_plus':
        #if date.day not in events:
            #continue
    # Generate names of daily directories
    dirday = os.path.join(dirAProf,date.strftime('%Y/%Y_%m_%d'))
    dirdayL1 = os.path.join(dirL1,date.strftime('%Y/%Y_%m_%d'))
    print(dirday)
    # List the content of the daily aeorosol directory
    if int(initials[-4:])>2019:
        ficL1 = sorted(glob.glob(dirdayL1+'/CAL_LID_L1-ValStage1-V3-40.*.hdf'))
    else:
        ficL1 = sorted(glob.glob(dirdayL1+'/CAL_LID_L1-Standard-V4-10.*.hdf'))
    print(len(ficL1))
    ficL1.reverse()
    # process all the half-orbits in the aerosol directory
    for i in range(len(ficL1)):
        # pop file
        file = ficL1.pop()
        print(file)
        # skip day files
        #if 'ZD' in file: continue
        # open file
        try:
            hdfL1 = SD(file)
        except :
            print('HDF4 Error -> giving up')
            continue
        # select orbits that intersect the box
        lonsL1 = hdfL1.select('Longitude').get()[:] % 360
        latsL1 = hdfL1.select('Latitude').get()[:]
        sel1L1 = (latsL1 < box[3]) & (latsL1 > box[2]) & (lonsL1 < box[1]) & (lonsL1 > box[0])
        # for anomalous cases where the orbit does not have a tropical segment
        if np.sum(sel1L1) == 0:
            print('orbit with no selected band')
            continue
        tai = hdfL1.select('Profile_Time').get()[:]
        tt = Time('1993-01-01 00:00:00',scale='tai') + TimeDelta(tai, format='sec')
        utc = tt.utc.datetime
        utcc = utc[sel1L1][0]+0.5*(utc[sel1L1][-1]-utc[sel1L1][0])
        fname = os.path.basename(file)
        # extract rootname
        rootname = fname[26:-4]
        #print(rootname,utcc,int(np.sum(sel1L1)),len(sel1L1))
        try:
            if int(initials[-4:])>2019:
                file = os.path.join(dirday,'CAL_LID_L2_05kmAPro-Prov-V3-40.'+rootname+'.hdf')
            else:
                file = os.path.join(dirday,'CAL_LID_L2_05kmAPro-Standard-V4-10.'+rootname+'.hdf')
            hdf = SD(file)
            # extract information form L1 file
            lats = hdf.select('Latitude').get()[:,1]
            lons = hdf.select('Longitude').get()[:,1] % 360
            sel1 = (lats < box[3]) & (lats > box[2]) & (lons < box[1]) & (lons > box[0])
            print('Apro',np.sum(sel1),len(sel1))
        except:
            print(file)
            print('Apro file cannot be found or opened')
            sel1 = None

        # generates a new record for this trace
        idx += 1
        print('idx',idx)
        Cald[idx] = {'fname':rootname,'utc':utcc,'sel1':sel1,'sel1L1':sel1L1}
        # store the dictionary of traces
        with gzip.open(namecal,'wb') as f:
            pickle.dump([box,date00,nbday+nbday0,Cald],f)

# # store the dictionary of traces
# with gzip.open(namecal,'wb') as f:
#     pickle.dump([box,date00,nbday+nbday0,Cald],f)

print('Now the Cald dictionary is stored in the directory :'+diroutput+' ; You can process to the etape2 program')
