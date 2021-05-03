#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Created on Sun May 31 21:45:54 2020

@author: Bernard Legras
"""
import numpy as np
from datetime import datetime,timedelta
from ECMWF_N import ECMWF

dat = ECMWF('FULL-EA',datetime(2020,1,21,6),exp=['F12','VOZ'])
dat._get_var('VO')
dat._get_var('U')
dat._get_var('V')
dat._get_var('T')
dat._get_var('O3')
# get forcast fields for same time
dat._get_var('VOF')
dat._get_var('UF')
dat._get_var('VF')
dat._get_var('TF')
dat._get_var('O3F')
# calculate PV for both
dat._mkp()
dat._mkthet()
dat._mkpv()
dat._mkz()
dat._mkp(suffix='F')
dat._mkthet(suffix='F')
dat._mkpv(suffix='F')

dat._mkinc(varss = ['PV','VO','O3','T'])

#%%
PVmean = np.mean(dat.var['PV'],axis=2)
O3mean = np.mean(dat.var['O3'],axis=2)
Tmean = np.mean(dat.var['T'],axis=2)
dat.var['PVano'] = dat.var['PV'] - PVmean[:,:,np.newaxis]
dat.var['Tano'] = dat.var['T'] - Tmean[:,:,np.newaxis]
dat.var['O3ano'] = dat.var['O3'] - O3mean[:,:,np.newaxis]
dats = dat.extract(varss=['Z','VO','VOD','PV','PVano','PVD','T','Tano','TD','O3','O3ano','O3D'],
                   lonRange=(240,300),latRange=(-80,-30))
dats.show('PV',48,figsize=(6,4))
dats.show('PVD',48,figsize=(6,4))
dats.show('O3',48,figsize=(6,4))
dats.show('O3D',48,figsize=(6,4))
dats.show('T',48,figsize=(6,4))
dats.show('TD',48,figsize=(6,4))
#%%
dats.chartlonz('PVano',-59,levs=(35,58),txt='PV anomaly (PVU)',scale=1.e6)
dats.chartlonz('PVD',-59,levs=(35,58),txt='PV increment (PVU)',scale=1.e6)
dats.chartlonz('Tano',-59,levs=(35,58),txt='T anomaly')
dats.chartlonz('TD',-59,levs=(35,58),txt='T increment')
dats.chartlonz('O3ano',-59,levs=(35,58),txt='O3 anomaly (mg/kg)',scale=1.e6)
dats.chartlonz('O3D',-59,levs=(35,58),txt='O3 increment (mg.kg)',scale=1.e6)
dats.chartlonz('VO',-59,levs=(35,58),txt=r'vorticity (10$^5$ s$^{-1})',scale=1.e5)
dats.chartlonz('VOD',-59,levs=(35,58),txt='vorticity increment (10$^5$ s$^{-1})',scale=1.e5)
