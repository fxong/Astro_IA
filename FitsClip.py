"""
Purpose:
    clip the fitsfile into sub-fitsfile with the desired area
Usage:
Caution:
    fitsfile may be fixed by CasaFix.py first
    the unit of x axis, y axis (and z axis) should refer to the head of fitsfile
History:
    Completed. Version 0.8.0. Oct. 08, 18.
    Axis infomation updated. Version 0.8.5. Oct. 18, 18.
    2D datacube support. TBC...
Copyright:
    written by fxong@CfA
"""

import os
pwd=os.getcwd()+'/'
os.chdir(pwd)
import numpy as np
import astropy.io.fits as ft
from astropy.wcs import WCS as wcs

fitsfile='' # without suffix
x=[] # deg or km/s
y=[] # deg or km/s
z=[] # deg or km/s
olshow=False # True/False

dim=[x,y,z]
hdr=ft.getheader(fitsfile+'.fits')
unit=[hdr['cunit1'], hdr['cunit2'], hdr['cunit3']]
for i in range(len(unit)):
    if unit[i]=='km/s' or unit[i]=='m/s':
        inax=i
        fac=1e3 if unit[i]=='km/s' else 1
    else:
        inax=False
if inax: dim[inax]=[i*1e3 for i in dim[inax]]

for i in range(len(dim)):
    if dim[i]!=[]:
        world=[hdr['crval1'], hdr['crval2'], hdr['crval3']]
        if inax: world[inax]=world[inax]*fac
        world[i]=dim[i]
        pixel=wcs(hdr).wcs_world2pix(world[0], world[1], world[2], 0)
        dim[i]=[int(round(pixel[i][0])), int(round(pixel[i][1]))]
        dim[i].sort()
    else:
        dim[i]=[0, hdr['naxis'+str(i+1)]-1]

if not olshow:
    dim.reverse()
    shape=[i[1]-i[0]+1 for i in dim]
    dlpix=[i[0] for i in dim]
    subdat=np.zeros(tuple(shape))
    dat=ft.getdata(fitsfile+'.fits')
    for i in range(shape[0]):
        for j in range(shape[1]):
            for k in range(shape[2]):
                subdat[i,j,k]=dat[i+dlpix[0],j+dlpix[1],k+dlpix[2]]

    rfval=wcs(hdr).wcs_pix2world(dlpix[2], dlpix[1], dlpix[0], 0)
    rfval=[float(i) for i in rfval]
    if inax: rfval[inax]=rfval[inax]/fac
    for i in range(len(rfval)):
        hdr['crval'+str(i+1)]=rfval[i]
        hdr['crpix'+str(i+1)]=float(1.0)

    if inax and fac==1:
        hdr['crval'+str(inax+1)]/=1e3
        hdr['cdelt'+str(inax+1)]/=1e3
        hdr['cunit'+str(inax+1)]='km/s'
    hdr['history']='Clipped from '+fitsfile+'.fits'

    ft.writeto(fitsfile+'_sub.fits', subdat, hdr, output_verify='fix+warn', overwrite=True)
else:
    print('The X range of the desired area is '+str(dim[0])+' in pixel.')
    print('The Y range of the desired area is '+str(dim[1])+' in pixel.')
    print('The Z range of the desired area is '+str(dim[2])+' in pixel.')
    print('The datacube starts with pixel 0.')
