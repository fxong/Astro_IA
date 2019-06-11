"""
Purpose:
    build and fit the radial density profiles of filaments
Usage:
Caution:
History:
    Completed. Version 0.9.0. Dec. 21, 18.
    Updated. Version 0.9.5. Feb. 20, 19.
Copyright:
    written by fxong@CfA
"""

import os
pwd=os.getcwd()+'/'
os.chdir(pwd)
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import ascii
from astropy.io import fits as ft
from astropy.table import Table
from astropy.wcs import WCS as wcs
from radfil import radfil_class, styles

fitsfile='' # column density fits file, without suffix
fits_err=None # error fits file of column density, without suffix, default is None
pathfile='' # path file of the filament, without suffix
distance=1 # distance to the source, unit is kpc, default is 1
dis_error=0 # error of the distance, unit is kpc, default is 0
step=3 # sample along the filament when cutting the profile, unit is pixel, default is 3
width=10 # width of the mask in which the peak pixels will be confined, unit is pixel, default is 10
shift=True # whether the center of the profile should be shifted to the peak pixel, defaut is True
binsz=None # size of the bin when averaging the profile, unit is pixel, default is 1.0; if no average is needed, set as None
xzoom=1.0 # factor to zoom in the x-axis of profile, set between 0.1 and 1.0, default is 1.0

fitfunc=None # fitting function, 'Plummer', 'Gaussian' or 'Both', default is None
bgdegree=0 # polynomial order used in background subtraction, 0 or 1, default is 0
bgdist=None # radial range used in background subtraction, like (0.3,0.5), unit is pc, default is None
fitdist_pm=0.1 # radial range in which to fit the profile using 'Plummer', unit is pc, default is 0.1
fitdist_gs=0.1 # radial range in which to fit the profile using 'Gaussian', unit is pc, default is 0.1

""" inputs """
fil_image, fil_header=ft.getdata(fitsfile+'.fits', header=True)
if fits_err is not None: err_image, err_header=ft.getdata(fits_err+'.fits', header=True)
fil_path=ascii.read(pathfile+'.path', format='fixed_width', names=['x', 'y'])
fil_path_x, fil_path_y=[i for i in fil_path['x']], [j for j in fil_path['y']]

path_pixel=wcs(fil_header).wcs_world2pix(fil_path_x, fil_path_y, 0)
fil_path=np.zeros_like(fil_image)
for i in range(path_pixel[0].shape[0]):
    fil_path[int(round(path_pixel[1][i])), int(round(path_pixel[0][i]))]=1.0
fil_path=fil_path.astype(bool)

major_fwhm=fil_header['bmaj']*3600
minor_fwhm=fil_header['bmin']*3600
beam_fwhm=np.sqrt(major_fwhm*minor_fwhm)

""" mask and sample """
radobjmk=radfil_class.radfil(fil_image, filspine=fil_path, header=fil_header, beamwidth=beam_fwhm, distance=distance*1e3)
radobjmk.build_profile(samp_int=step, cutdist=width*radobjmk.imgscale.value)
plt.close()
if binsz is not None:
    bins=int(round(abs(np.nanmax(radobjmk.xall)-np.nanmin(radobjmk.xall))/(radobjmk.imgscale.value*binsz)))
else: bins=None

radobj=radfil_class.radfil(fil_image, mask=radobjmk.mask, filspine=fil_path, header=fil_header, beamwidth=beam_fwhm, distance=distance*1e3)
radobj.build_profile(samp_int=step, bins=bins, shift=shift)
plt.savefig('FilProfile_cuttings.png', format='png')
if fitfunc is None: plt.show()
else: plt.close()

if fits_err is not None:
    radobj_err=radfil_class.radfil(err_image, mask=radobjmk.mask, filspine=fil_path, header=err_header, beamwidth=beam_fwhm, distance=distance*1e3)
    radobj_err.build_profile(samp_int=step, bins=bins, shift=shift)
    plt.close()

fil_mask=radobj.mask.astype('float32')
index=np.where(fil_path==True)
for i in range(index[0].shape[0]):
    fil_mask[index[0][i], index[1][i]]=2.0
ft.writeto('FilProfile_mask.fits', fil_mask, fil_header, output_verify='fix+warn', overwrite=True)
for i in range(index[0].shape[0]):
    fil_mask[index[0][i], index[1][i]]=np.nan
ft.writeto('FilProfile_mask_image.fits', fil_mask*fil_image, fil_header, output_verify='fix+warn', overwrite=True)

if os.path.isfile('./FilProfile_fitting.data'): os.system('rm -rf FilProfile_fitting.data')
info=open('FilProfile_fitting.data', mode='w')
fil_length=index[0].shape[0]*radobj.imgscale.value
fil_length_error=(dis_error*1e3)*(fil_length/(distance*1e3))
print('\nThe length of the filament is %f +/- %f pc' %(fil_length, fil_length_error))
info.write('\nThe length of the filament is %f +/- %f pc\n' %(fil_length, fil_length_error))
mask_width=np.nanmedian(radobj.dictionary_cuts['mask_width'])
print('The average width of the mask is %f pc\n' %mask_width)
info.write('The average width of the mask is %f pc\n' %mask_width)

""" built the profile """
if bins is not None:
    bins=np.linspace(np.nanmin(radobj.xall), np.nanmax(radobj.xall), bins+1)
    binsy=[radobj.yall[((radobj.xall >= (i-0.5*np.diff(bins)[0])) & (radobj.xall < (i+0.5*np.diff(bins)[0])))] for i in radobj.masterx]
    stdevy=np.array([np.nanstd(i) for i in binsy])
    if fits_err is not None:
        binsy_err=[radobj_err.yall[((radobj_err.xall >= (i-0.5*np.diff(bins)[0])) & (radobj_err.xall < (i+0.5*np.diff(bins)[0])))] for i in radobj_err.masterx]
        erry=np.array([np.sqrt(np.nanmean(np.power(i,2))) for i in binsy_err])
    plt.errorbar(radobj.masterx, radobj.mastery, color='black', yerr=stdevy, ecolor='gold', capsize=2)
else: plt.plot(radobj.masterx, radobj.mastery, color='black', alpha=0.6)

plt.xlabel('Radial Distance (pc)', size=20)
plt.ylabel(fil_header['btype']+' ('+fil_header['bunit']+')', size=20)
plt.tick_params(labelsize=20)
xlim=[abs(np.nanmin(radobj.xall)), abs(np.nanmax(radobj.xall))]
xlim=xzoom*np.nanmin(xlim)
plt.xlim(-xlim, xlim)
plt.plot(plt.xlim(), [0,0], linestyle='dashed')

if bins is not None:
    xran=((radobj.masterx >= -xlim) & (radobj.masterx <= xlim))
    text='bin size = %.1f pixel(s)\nnumber of bins = %i\npoints in each bin: %i~%i' %(binsz, np.sum(xran), np.nanmin(radobj.masternobs[xran]), np.nanmax(radobj.masternobs[xran]))
    plt.text(0.03, 0.95, text, weight='bold', size=20, verticalalignment='top', transform=plt.gca().transAxes)
else:
    text='number of profiles = %i' %len(radobj.dictionary_cuts['profile'])
    plt.text(0.03, 0.95, text, weight='bold', size=20, verticalalignment='top', transform=plt.gca().transAxes)
plt.savefig('FilProfile_profile.png', format='png')
if fitfunc is None: plt.show()
else: plt.close()

""" Plummer fitting """
if fitfunc=='Plummer' or fitfunc=='Both':
    radobj.fit_profile(fitfunc='Plummer', fitdist=fitdist_pm, bgdegree=bgdegree, bgdist=bgdist, beamwidth=beam_fwhm, verbose=False)

    if bgdist is not None:
        if bins is not None:
            plt.errorbar(radobj.masterx, radobj.mastery, fmt='none', yerr=stdevy, ecolor='gold', capsize=2)
        else: plt.plot(radobj.masterx, radobj.mastery, color='black', alpha=0.6)
        radobj.plotter().plotFits(plt.gca(), 'bg')
        plt.xlabel('Radial Distance (pc)', size=20)
        plt.ylabel(fil_header['btype']+' ('+fil_header['bunit']+')', size=20)
        plt.tick_params(labelsize=20)
        plt.xlim(-xlim, xlim)
        plt.plot(plt.xlim(), [0,0], linestyle='dashed')
        plt.savefig('FilProfile_Plummer_background.png', format='png')
        plt.show()
        ploty=radobj.mastery-radobj.bgfit(radobj.masterx)
    else: ploty=radobj.mastery

    if bins is not None:
        plt.errorbar(radobj.masterx, ploty, fmt='none', yerr=stdevy, ecolor='gold', capsize=2)
    else: plt.plot(radobj.masterx, ploty, color='black', alpha=0.6)
    radobj.plotter().plotFits(plt.gca(), 'model')
    plt.xlabel('Radial Distance (pc)', size=20)
    plt.ylabel(fil_header['btype']+' ('+fil_header['bunit']+')', size=20)
    plt.tick_params(labelsize=20)
    plt.xlim(-xlim, xlim)
    plt.plot(plt.xlim(), [0,0], linestyle='dashed')
    plt.savefig('FilProfile_Plummer_fitting.png', format='png')
    plt.show()

    print('\n==== Plummer-like Fitting ====')
    info.write('\n==== Plummer-like Fitting ====\n')
    for (name, value, error) in zip(radobj.profilefit.param_names, radobj.profilefit.parameters, radobj.std_error):
        print('The best-fit %s is %f +/- %f' %(name, value, error))
        info.write('The best-fit %s is %f +/- %f\n' %(name, value, error))
    # print('The non-deconvolved FWHM is %f pc' %radobj.FWHM)
    # info.write('The non-deconvolved FWHM is %f pc\n' %radobj.FWHM)
    # beamwidth_phys=(radobj.beamwidth/radobj.imgscale_ang).decompose()*radobj.imgscale.value
    # print('The physical size of the beam is %f pc' %beamwidth_phys)
    # info.write('The physical size of the beam is %f pc\n' %beamwidth_phys)
    # print('The deconvolved FWHM is %f pc\n' %radobj.FWHM_deconv)
    # info.write('The deconvolved FWHM is %f pc\n' %radobj.FWHM_deconv)

""" Gaussian fitting """
if fitfunc=='Gaussian' or fitfunc=='Both':
    radobj.fit_profile(fitfunc='Gaussian', fitdist=fitdist_gs, bgdegree=bgdegree, bgdist=bgdist, beamwidth=beam_fwhm, fix_mean=True, verbose=False)

    if bgdist is not None:
        if bins is not None:
            plt.errorbar(radobj.masterx, radobj.mastery, fmt='none', yerr=stdevy, ecolor='gold', capsize=2)
        else: plt.plot(radobj.masterx, radobj.mastery, color='black', alpha=0.6)
        radobj.plotter().plotFits(plt.gca(), 'bg')
        plt.xlabel('Radial Distance (pc)', size=20)
        plt.ylabel(fil_header['btype']+' ('+fil_header['bunit']+')', size=20)
        plt.tick_params(labelsize=20)
        plt.xlim(-xlim, xlim)
        plt.plot(plt.xlim(), [0,0], linestyle='dashed')
        plt.savefig('FilProfile_Gaussian_background.png', format='png')
        plt.show()
        ploty=radobj.mastery-radobj.bgfit(radobj.masterx)
    else: ploty=radobj.mastery

    if bins is not None:
        plt.errorbar(radobj.masterx, ploty, fmt='none', yerr=stdevy, ecolor='gold', capsize=2)
    else: plt.plot(radobj.masterx, ploty, color='black', alpha=0.6)
    radobj.plotter().plotFits(plt.gca(), 'model')
    plt.xlabel('Radial Distance (pc)', size=20)
    plt.ylabel(fil_header['btype']+' ('+fil_header['bunit']+')', size=20)
    plt.tick_params(labelsize=20)
    plt.xlim(-xlim, xlim)
    plt.plot(plt.xlim(), [0,0], linestyle='dashed')
    plt.savefig('FilProfile_Gaussian_fitting.png', format='png')
    plt.show()

    print('\n==== Gaussian Fitting ====')
    info.write('\n==== Gaussian Fitting ====\n')
    if radobj.std_error.shape[0]==2: radobj.std_error=np.insert(radobj.std_error, 1, 0.0)
    for (name, value, error) in zip(radobj.profilefit.param_names, radobj.profilefit.parameters, radobj.std_error):
        print('The best-fit %s is %f +/- %f' %(name, value, error))
        info.write('The best-fit %s is %f +/- %f\n' %(name, value, error))
    gau_std=radobj.profilefit.parameters[2]
    gau_std_err=radobj.std_error[2]
    gau_FWHM=2*np.sqrt(2*np.log(2))*gau_std
    gau_FWHM_err=2*np.sqrt(2*np.log(2))*gau_std_err
    print('The non-deconvolved FWHM is %f +/- %f pc' %(radobj.FWHM, gau_FWHM_err))
    info.write('The non-deconvolved FWHM is %f +/- %f pc\n' %(radobj.FWHM, gau_FWHM_err))
    beamwidth_phys=(radobj.beamwidth/radobj.imgscale_ang).decompose()*radobj.imgscale.value
    print('The physical size of the beam is %f pc' %beamwidth_phys)
    info.write('The physical size of the beam is %f pc\n' %beamwidth_phys)
    gau_FWHM_dec=np.sqrt(np.power(gau_FWHM,2)-np.power(beamwidth_phys,2))
    gau_FWHM_dec_err=gau_FWHM*gau_FWHM_err/gau_FWHM_dec
    print('The deconvolved FWHM is %f +/- %f pc\n' %(radobj.FWHM_deconv, gau_FWHM_dec_err))
    info.write('The deconvolved FWHM is %f +/- %f pc\n' %(radobj.FWHM_deconv, gau_FWHM_dec_err))

info.close()

""" outputs """
if (bins is not None) and (fitfunc is not None):
    if fits_err is None:
        erry=np.array([i for i in stdevy])
        erry[np.invert(np.isnan(erry))]=float(0)
    tabxy=Table([radobj.masterx, ploty, stdevy, erry], names=['Radial Distance (pc)', fil_header['btype']+' ('+fil_header['bunit']+')', 'Standard Deviation'+' ('+fil_header['bunit']+')', 'Uncertainties'+' ('+fil_header['bunit']+')'])
    if bgdist is None: output='FilProfile_profile.data'
    else: output='FilProfile_profile_bgrm.data'
    ascii.write(tabxy, output=output, format='fixed_width', overwrite=True)

if os.path.isdir('./FilProfile'): os.system('rm -rf FilProfile')
os.mkdir('FilProfile')
os.system('mv FilProfile*.fits FilProfile/')
os.system('mv FilProfile*.png FilProfile/')
os.system('mv FilProfile*.data FilProfile/')
