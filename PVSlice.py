"""
Purpose:
    extract the position-velocity slice from a fitsfile
Usage:
Cautious:
    ds9 region should be in Line shape
History:
    Completed. Version 0.5.0. Oct. 18, 18.
    Path export. TBC...
Copyright: written by fxong@CfA
"""

import os
pwd=os.getcwd()+'/'
os.chdir(pwd)
from astropy import units as u
from astropy.coordinates import SkyCoord
from regions import read_ds9
from spectral_cube import SpectralCube
from pvextractor import Path
from pvextractor import extract_pv_slice

fitsfile='' # without suffix
spec=[] # in km/s
ds9path='' # ds9 region file without suffix
pathx=[] # in deg
pathy=[] # in deg
step=1 # in pixel, defaut is 1
width=0 # in arcsec, defaut is 0
frame='icrs' # icrs, galactic

""" cube and path """
cube=SpectralCube.read(fitsfile+'.fits')
if spec!=[]:
    cube=cube.spectral_slab(spec[0]*u.km/u.s, spec[1]*u.km/u.s)

if ds9path!='':
    regs=read_ds9(ds9path+'.reg')
    xy=list(range(len(regs)*2))
    for i in range(len(regs)):
        if frame=='icrs':
            xy[2*i]=[regs[i].start.ra.value, regs[i].start.dec.value]
            xy[2*i+1]=[regs[i].end.ra.value, regs[i].end.dec.value]
        elif frame=='galactic':
            xy[2*i]=[regs[i].start.l.value, regs[i].start.b.value]
            xy[2*i+1]=[regs[i].end.l.value, regs[i].end.b.value]
    xy.sort()
    pathx, pathy=list(range(len(xy))), list(range(len(xy)))
    for i in range(len(xy)):
        pathx[i], pathy[i]=xy[i][0], xy[i][1]

value=SkyCoord(pathx, pathy, frame=frame, unit=u.deg)
if width!=0:
    path=Path(value, width=width*u.arcsec)
else:
    path=Path(value)

""" path resample """
linex, liney=path.sample_points(spacing=step, wcs=cube.wcs)
linexy=cube.wcs.wcs_pix2world(linex, liney, 0, 0)
linex=[i for i in linexy[0]]
liney=[i for i in linexy[1]]

if width!=0:
    polyxy=path.sample_polygons(spacing=step, wcs=cube.wcs)
    polyx, polyy=[], []
    for i in range(len(polyxy)):
        polyx+=polyxy[i].x
        polyy+=polyxy[i].y
    polyxy=cube.wcs.wcs_pix2world(polyx, polyy, 0, 0)
    polyx=[i for i in polyxy[0]]
    polyy=[i for i in polyxy[1]]

""" slice export """
slice=extract_pv_slice(cube, path, spacing=step)

slice.data=slice.data.transpose()
keys=['CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CUNIT']
for i in range(len(keys)):
    slice.header[keys[i]+'1'], slice.header[keys[i]+'2']=slice.header[keys[i]+'2'], slice.header[keys[i]+'1']
    slice.header.comments[keys[i]+'1'], slice.header.comments[keys[i]+'2']=slice.header.comments[keys[i]+'2'], slice.header.comments[keys[i]+'1']

slice.header['CDELT1']=slice.header['CDELT1']/1e3
slice.header.comments['CDELT1']='[km/s] Coordinate increment at reference point'
slice.header['CUNIT1']='km/s'
slice.header['CRVAL1']=slice.header['CRVAL1']/1e3
slice.header.comments['CRVAL1']='[km/s] Coordinate value at reference point'
if spec!=[]:
    slice.header['CRPIX1']=float(1.0)
    slice.header['CRVAL1']=linexy[2][0]/1e3

slice.writeto(fitsfile+'_'+ds9path+'.fits',output_verify='fix+warn', overwrite=True)

""" path info export """
# from astropy.io import ascii
# from astropy.table import Table
