"""
Purpose:
    fix the fitsfile exported by CASA, remove the redundant axis and change the unit of velocity axis
Usage:
Cautious:
History:
    Completed. Version 0.1.0. Sep. 14, 18.
    Updated. Version 0.8.0. Sep. 20, 18.
Copyright:
    written by fxong@CfA
"""

import os
pwd=os.getcwd()+'/'
os.chdir(pwd)
import astropy.io.fits as ft
from numpy import reshape as rs
# import astropy.units as u

fitsfile='' # without suffix
hdul=ft.open(fitsfile+'.fits')
hdr, dat=hdul[0].header, hdul[0].data
hdul.close()

""" remove redundant axis """
if hdr['NAXIS']==4 and hdr['CTYPE4']=='STOKES':
    oldsp=list(dat.shape)
    oldsp.pop(hdr['NAXIS']-4)
    newsp=tuple(oldsp)
    dat=rs(dat,newsp)
    
    del hdr['CTYPE4'], hdr['CRVAL4'], hdr['CDELT4'], hdr['CRPIX4'], hdr['CUNIT4']
    del hdr['PC4_1'], hdr['PC4_2'], hdr['PC4_3'], hdr['PC1_4'], hdr['PC2_4'], hdr['PC3_4'], hdr['PC4_4']

""" change velocity unit """
if hdr['CUNIT3']=='m/s':
    hdr['CRVAL3']=hdr['CRVAL3']/1e3
    hdr['CDELT3']=hdr['CDELT3']/1e3
    hdr['CUNIT3']='km/s'

hdr['HISTORY']='Modified from '+fitsfile+'.fits'

# fitsfile='' # without suffix
ft.writeto(fitsfile+'_rdv.fits', dat, hdr, output_verify='fix+warn', overwrite=True)