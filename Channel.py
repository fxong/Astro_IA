"""
Purpose:
    extract the position-velocity slice from a fitsfile
Usage:
Caution:
    fitsfile may be fixed by CasaFix.py first
    the unit of 'spec' is km/s
History:
    Completed. Version 0.8.0. Oct. 26, 18.
Copyright:
    written by fxong@CfA
"""

import os
pwd=os.getcwd()+'/H13COP/FoF'
os.chdir(pwd)
import numpy as np
from astropy.io import fits as ft
from spectral_cube import SpectralCube as sc
from astropy.wcs import WCS as wcs

fitsfile='H13COP_pbcor_rdv_sub' # without suffix
spec=[-12.318,2.081] # in km/s
nchan=20 # in channel
rms=0.070 # in Jy/beam, default is 0
nexc=2 # pixels with value < -next*rms are excluded when doing moment0, default is 3/1.5
ncut=3 # pixel with value < ncut*integrated rms are set to be zero after doing moment0, default is 3/5

""" input normalization """
spec=[i*1e3 for i in spec]
dat, hdr=ft.getdata(fitsfile+'.fits', header=True)
pixel=wcs(hdr).wcs_world2pix(hdr['crval1'], hdr['crval2'], spec, 0)
chans=[round(pixel[2][0]), round(pixel[2][1])]
dchan=round(abs(chans[1]-chans[0])/nchan)
for i in range(nchan):
    chans.insert(i+1, chans[i]+dchan)
chans.pop()

""" print information """
world=wcs(hdr).wcs_pix2world(hdr['crpix1'], hdr['crpix2'], chans, 0)
spec=[i/1e3 for i in world[2]]
dvelo=dchan*hdr['cdelt3']
print('Extraction of %i channels from %f km/s to %f km/s.' %(nchan, spec[0], spec[-1]))
print('The width of each channel is %f km/s (%i times of spectral resolution).' %(dvelo, dchan))
print('The ranges of each channel are:')
for i in range(len(spec)-1):
    print('%f km/s to %f km/s' %(spec[i], spec[i+1]))

""" moment0 of each channel """
for i in range(dat.shape[0]):
    for j in range(dat.shape[1]):
        for k in range(dat.shape[2]):
            if rms>0 and np.invert(np.isnan(dat[i,j,k])) and dat[i,j,k]<=-nexc*rms: dat[i,j,k]=float(0.0)
hdu=ft.PrimaryHDU(data=dat, header=hdr)
hdul=ft.HDUList([hdu])
cube=sc.read(hdul)
hdul.close()

shape=(nchan, dat.shape[1], dat.shape[2])
newdat=np.zeros(shape)
chans=[int(i) for i in chans]
for i in range(shape[0]):
    subcube=cube[chans[i]:chans[i+1],:,:]
    subcube=subcube.moment(order=0)
    newdat[i]=subcube.value
    if rms>0: Irms=rms*np.sqrt(dvelo*(dvelo/dchan))
    for j in range(shape[1]):
        for k in range(shape[2]):
            if rms>0 and newdat[i][j,k]<=ncut*Irms: newdat[i][j,k]=float(0.0)

""" output fitsfile """
hdr['btype']='Moment0'
hdr['bunit']='Jy/beam.km/s'
hdr['crval3']=np.mean([spec[0],spec[1]])
hdr['cdelt3']=dvelo
hdr['crpix3']=float(1.0)
hdr['cunit3']='km/s'
hdr['comment']='Extraction of %i channels from %f km/s to %f km/s.' %(nchan, spec[0], spec[-1])
hdr['comment']='The width of each channel is %f km/s (%i times of spectral resolution).' %(dvelo, dchan)
hdr['comment']='The ranges of each channel are:'
for i in range(len(spec)-1):
    hdr['comment']='%f km/s to %f km/s' %(spec[i], spec[i+1])

if rms>0:
    ft.writeto(fitsfile+'_channels_3rms.fits', newdat, hdr, output_verify='fix+warn', overwrite=True)
else:
    ft.writeto(fitsfile+'_channels.fits', newdat, hdr, output_verify='fix+warn', overwrite=True)