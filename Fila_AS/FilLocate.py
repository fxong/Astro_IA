"""
Purpose:
    locate the 2D positions of filaments within molecular clouds
Usage:
Caution:
    fitsfile should be cut off by rms threshold first
History:
    Completed. Version 0.5.0. Oct. 25, 18.
Copyright:
    written by fxong@CfA
"""

import os
pwd=os.getcwd()+'/'
os.chdir(pwd)

import numpy as np
import astropy.units as u
from astropy.io import fits as ft
from fil_finder import FilFinder2D
import matplotlib.pyplot as plt

fitsfile='' # without suffix
distance= # in kpc
scale=100 # default is 100 in %
smooth=0.05 # default is 0.05 in pc
width=0.1 # default is 0.1 in pc
aspect=5 # default is 5
branch=3 # default is 3 in pixel
save_dir='FilFinder' # output directory

""" inputs """
hdul=ft.open(fitsfile+'.fits')
hdr, dat=hdul[0].header, hdul[0].data
hdul.close()
fil=FilFinder2D(dat, header=hdr, distance=distance*u.kpc, save_name=save_dir)

""" rescale the image """
fil.preprocess_image(flatten_percent=scale)
plt.subplot(121)
plt.imshow(fil.image.value, origin='lower', cmap='binary')
plt.subplot(122)
plt.imshow(fil.flat_img.value, origin='lower', cmap='binary')
plt.savefig(save_dir+'_flat.png')
plt.show()

""" create mask and skeletons """
fil.create_mask(smooth_size=smooth*u.pc, adapt_thresh=width*u.pc, size_thresh=np.pi*aspect*(width*u.pc)**2, verbose=True, save_png=True)
fil.medskel(verbose=True, save_png=True)

""" prune the branches """
plt.imshow(fil.flat_img.value, origin='lower', cmap='binary')
plt.contour(fil.skeleton, colors='red')
fil.analyze_skeletons(prune_criteria='length', branch_thresh=branch*u.pix)
plt.contour(fil.skeleton, colors='lime')
plt.savefig(save_dir+'_final_skeletons.png')
plt.show()

""" outputs """
skedat=fil.skeleton
for i in range(skedat.shape[0]):
    for j in range(skedat.shape[1]):
        if np.isnan(dat[i,j]): skedat[i,j]=np.nan
ft.writeto(save_dir+'_skeletons.fits', skedat, hdr, output_verify='fix+warn', overwrite=True)

skedat=fil.skeleton
for i in range(dat.shape[0]):
    for j in range(dat.shape[1]):
        if np.invert(np.isnan(skedat[i,j])) and skedat[i,j]==1:
            dat[i,j]=np.nan
ft.writeto(save_dir+'_skeletons_image.fits', dat, hdr, output_verify='fix+warn', overwrite=True)

if os.path.isdir('./'+save_dir): os.system('rm -rf '+save_dir)
os.mkdir(save_dir)
os.system('mv '+save_dir+'*.png '+save_dir+'/')
