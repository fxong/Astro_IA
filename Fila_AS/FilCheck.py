"""
Purpose:
    verify the velocity coherence of filaments channel by channel within the data cube
Usage:
Caution:
History:
    Completed. Version 0.5.0. Nov. 25, 18.
Copyright:
    written by fxong@CfA
"""

import os
pwd=os.getcwd()+'/'
os.chdir(pwd)

import numpy as np
from astropy.io import fits as ft
from astropy.wcs import WCS as wcs

fitsfile='' # without suffix
filafile='' # without suffix
chans=[] # range of channels used in the check
rms=0 # Jy/beam, default is 0

dat, hdr=ft.getdata(fitsfile+'.fits', header=True)
fdat, fhdr=ft.getdata(filafile+'.fits', header=True)

sdat=np.zeros_like(dat)
for i in range(chans[0],chans[1]+1):
    for j in range(dat.shape[1]):
        for k in range(dat.shape[2]):
            if fdat[j,k]==1:
                if dat[i,j,k]>3*rms:
                    dat[i,j,k]=np.nan
                    sdat[i,j,k]=float(1.0)

ft.writeto('FilCheck_skeletons_image.fits', dat, hdr, output_verify='fix+warn', overwrite=True)
ft.writeto('FilCheck_skeletons.fits', sdat, hdr, output_verify='fix+warn', overwrite=True)
