"""
Purpose:
    generate the moment maps of fitsfile using CASA package
Usage:
    casa -c Moments.py
Cautious:
    arguments 'chans' and 'box' start with pixel 0
    Astropy package is needed to be installed in CASA
    fitsfile may be fixed by CasaFix.py first
    this script only works for M0(ra, dec, spec), M1(spec), M2(spec), M6(spec)
History:
    Completed. Version 0.1.0. Sep. 14, 18.
    Updated. Version 0.8.0. Sep. 20, 18.
    Algorithm changed. Version 0.9.0. Oct. 08, 18.
Copyright:
    written by fxong@CfA
"""

import os
import numpy as np
from math import pi
import astropy.io.fits as ft
from astropy.wcs import WCS as wcs

fitsfile='' # without suffix
moment=[] # M0, M1, M2, M6(rms)
axis='' # ra, dec, spec
chans=[] # [spec0,spec1(, spec2,spec3)] in channel
box=[] # [ra0,dec0,ra1,dec1(, ra2,dec2,ra3,dec3)] in pixel
rms=0 # Jy/beam
inpix=[-1] # [-1], [3*rms,1e3]
expix=[-1] # [-1], [-1e3,-3*rms/-1]
Nrms=0 # n*rms cutoff, default is 0

""" input normalization """
ohdr=ft.getheader(fitsfile+'.fits')

if chans!=[]:
    spec=chans
    numc, strc=(len(chans)/2), ''
    for i in range(int(numc)):
        strc=strc+str(chans[0+2*i])+'~'+str(chans[1+2*i])+';'
    chans=strc.rstrip(';')
else:
    spec=[0, ohdr['NAXIS3']-1]
    chans=''

if box!=[]:
    numra, ra=(len(box)/4)*2, []
    for i in range(int(numra)): ra.append(box[0+2*i])
    numdec, dec=(len(box)/4)*2, []
    for i in range(int(numdec)): dec.append(box[1+2*i])
    box=[str(i) for i in box]
    box=','.join(box)
else:
    ra=[0, ohdr['NAXIS1']-1]
    dec=[0, ohdr['NAXIS2']-1]
    box=''

""" integrated range retrieval """
if axis=='ra': inax=1
elif axis=='dec': inax=2
elif axis=='spec': inax=3

dim=[ra,dec,spec]
pixel=[ohdr['CRPIX1'], ohdr['CRPIX2'], ohdr['CRPIX3']]
pixel[inax-1]=dim[inax-1]
rlpix=pixel[inax-1]
rlval=wcs(ohdr).wcs_pix2world(pixel[0], pixel[1], pixel[2], 0)
rlval[2]=rlval[2]/1e3
rlval=list(rlval[inax-1])

""" moment """
importfits(fitsfile+'.fits', fitsfile+'.im', zeroblanks=True, overwrite=True)

if os.path.isdir('./'+fitsfile+'_mom.im'):
    os.system('rm -f '+fitsfile+'_mom.im')
immoments(fitsfile+'.im', moments=moment, axis=axis, chans=chans, box=box, includepix=inpix, excludepix=expix, outfile=fitsfile+'_mom.im')

""" export and header fix """
exportfits(fitsfile+'_mom.im', fitsfile+'_mom.fits', velocity=True, dropstokes=True, stokeslast=False, overwrite=True)
tdat, thdr=ft.getdata(fitsfile+'_mom.fits', header=True)

if thdr['CUNIT3']=='m/s':
    thdr['CRVAL3']=thdr['CRVAL3']/1e3
    thdr['CDELT3']=thdr['CDELT3']/1e3
    thdr['CUNIT3']='km/s'

if thdr['NAXIS'+str(inax)]==1:
    oldsp=list(tdat.shape)
    oldsp.pop(thdr['NAXIS']-inax)
    newsp=tuple(oldsp)
    tdat=np.reshape(tdat,newsp)
    if axis=='ra': tdat=tdat.transpose()

    # if axis=='ra':
    #     thdr['CTYPE1'], thdr['CRVAL1'], thdr['CDELT1'], thdr['CRPIX1'], thdr['CUNIT1']=thdr['CTYPE2'], thdr['CRVAL2'], thdr['CDELT2'], thdr['CRPIX2'], thdr['CUNIT2']
    #     thdr['CTYPE2'], thdr['CRVAL2'], thdr['CDELT2'], thdr['CRPIX2'], thdr['CUNIT2']=thdr['CTYPE3'], thdr['CRVAL3'], thdr['CDELT3'], thdr['CRPIX3'], thdr['CUNIT3']
    # elif axis=='dec':
    #     thdr['CTYPE2'], thdr['CRVAL2'], thdr['CDELT2'], thdr['CRPIX2'], thdr['CUNIT2']=thdr['CTYPE3'], thdr['CRVAL3'], thdr['CDELT3'], thdr['CRPIX3'], thdr['CUNIT3']
    # else:
    #     pass

    thdr['CTYPE'+str(inax)], thdr['CRVAL'+str(inax)], thdr['CDELT'+str(inax)], thdr['CRPIX'+str(inax)], thdr['CUNIT'+str(inax)]=thdr['CTYPE3'], thdr['CRVAL3'], thdr['CDELT3'], thdr['CRPIX3'], thdr['CUNIT3']

    del thdr['CTYPE3'], thdr['CRVAL3'], thdr['CDELT3'], thdr['CRPIX3'], thdr['CUNIT3']
    del thdr['PC3_1'], thdr['PC3_2'], thdr['PC1_3'], thdr['PC2_3'], thdr['PC3_3']

thdr['BTYPE']='Moment'+str(moment[0])
unit=' km/s' if axis=='spec' else ' deg'
numval, strval=(len(rlval)/2), ''
for i in range(int(numval)):
    strval=strval+str(rlval[0+2*i])+'~'+str(rlval[1+2*i])+'; '
thdr['COMMENT']='Moment'+str(moment[0])+' (Axis: '+axis+') Integrated Range: '+'['+strval.rstrip('; ')+']'+unit

if thdr['BUNIT']=='Jy/beam.rad':
    tdat=tdat*(180./pi)
    thdr['BUNIT']='Jy/beam.deg'

""" n*rms cutoff """
if moment==[0] and Nrms>0 and axis=='spec':
    numf1, fac1=(len(rlval)/2), 0
    for i in range(int(numf1)):
        fac1=fac1+abs(rlval[1+2*i]-rlval[0+2*i])
    numf2, fac2=(len(rlpix)/2), 0
    for i in range(int(numf2)):
        fac2=fac2+abs(rlpix[1+2*i]-rlpix[0+2*i])
    Irms=rms*np.sqrt(fac1*(float(fac1)/float(fac2)))

    for i in range(tdat.shape[0]):
        for j in range(tdat.shape[1]):
            if tdat[i,j]<=Nrms*Irms: tdat[i,j]=0
    ft.writeto(fitsfile+'_m'+str(moment[0])+'_'+axis+'_3rms.fits', tdat, thdr, output_verify='fix+warn', overwrite=True)
    print('The pixels < %i*rms (%f %s) are set to be Zero' %(Nrms, Nrms*Irms, thdr['BUNIT']))
else:
    ft.writeto(fitsfile+'_m'+str(moment[0])+'_'+axis+'.fits', tdat, thdr, output_verify='fix+warn', overwrite=True)

os.system('rm -f '+fitsfile+'_mom.fits')
os.system('rm -rf *.im *.last *.log')
