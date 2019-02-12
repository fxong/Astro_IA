""" 
Purpose:
    calculate the optical depths of two isotopologues assuming they have the same excitation temperature
Usage:
Caution:
    the input fitsfiles must be 2D and should have already been applied by 3*rms cutoff
    the tracer in fitsfile1 is assumed to be optically thicker than the tracer in fitsfile2
    the optical depth of tracer in fitsfile2 should be in the range 0.00~2.00
History:
    Completed. Version 0.8.0. Nov. 07, 18.
    Updated. Version 0.9.0. Jan. 11, 19.
Copyright:
    written by fxong@CfA
"""

from astropy.io import fits as ft
import numpy as np

fitsfile1='' # without suffix
fitsfile2='' # without suffix
ratio='CO' # optical depth ratio, [tracer in fitsfile1]/[tracer in fitsfile2]; input a number or the string 'CO' (only available for 12CO/13CO)
distance=1.0 # kpc, distance is needed when ratio equals 'CO'
step=0.001

dat1, hdr1=ft.getdata(fitsfile1+'.fits', header=True)
dat2, hdr2=ft.getdata(fitsfile2+'.fits', header=True)

if ratio=='CO':
    Dsun=8.34 # Reid et al.(2014)
    Ratio_12CO_13CO=5.41*(Dsun+distance)+19.03 # Milam et al.(2005)
    ratio=Ratio_12CO_13CO

for i in range(dat1.shape[0]):
    for j in range(dat1.shape[1]):
        if dat1[i,j]==0: dat1[i,j]=np.nan
for i in range(dat2.shape[0]):
    for j in range(dat2.shape[1]):
        if dat2[i,j]==0: dat2[i,j]=np.nan
dat1_dat2=dat1/dat2
minimum=np.nanmin(dat1_dat2)
maximum=np.nanmax(dat1_dat2)
mean=np.nanmean(dat1_dat2)
median=np.nanmedian(dat1_dat2)

tau2=np.arange(step, 2+step, step)
tau1=tau2*ratio
T1_T2=(1-np.exp(-tau1))/(1-np.exp(-tau2))

tau2_dat=np.zeros_like(dat1_dat2)
tau1_dat=np.zeros_like(dat1_dat2)
for i in range((dat1_dat2).shape[0]):
    for j in range((dat1_dat2).shape[1]):
        if np.isnan(dat1_dat2[i,j]):
            tau2_dat[i,j]=np.nan
            tau1_dat[i,j]=np.nan
        if np.invert(np.isnan(dat1_dat2[i,j])):
            index=np.argmin(abs(T1_T2-dat1_dat2[i,j]))
            tau2_dat[i,j]=tau2[index]
            tau1_dat[i,j]=tau1[index]
index1=np.argmin(abs(T1_T2-minimum))
index2=np.argmin(abs(T1_T2-maximum))
index3=np.argmin(abs(T1_T2-mean))
index4=np.argmin(abs(T1_T2-median))

print('brigthness temperature ratio (tracer in %s.fits / tracer in %s.fits):' %(fitsfile1, fitsfile2))
print('(1)min: %f; (2)max: %f; (3)mean: %f; (4)median: %f' %(minimum, maximum, mean, median))
print('\nWhen applying the min/max/mean/median brigthness temperature ratio, the corresponding optical depth of tracer in %s.fits:' %(fitsfile1))
print('(1): %f; (2): %f; (3): %f; (4): %f' %(tau1[index1], tau1[index2], tau1[index3], tau1[index4]))
print('\nWhen applying the min/max/mean/median brigthness temperature ratio, the corresponding optical depth of tracer in %s.fits:' %(fitsfile2))
print('(1): %f; (2): %f; (3): %f; (4): %f' %(tau2[index1], tau2[index2], tau2[index3], tau2[index4]))
print('\noptical depth of tracer in %s.fits:' %(fitsfile1))
print('min: %f; max: %f; mean: %f; median: %f' %(np.nanmin(tau1_dat), np.nanmax(tau1_dat), np.nanmean(tau1_dat), np.nanmedian(tau1_dat)))
print('\noptical depth of tracer in %s.fits:' %(fitsfile2))
print('min: %f; max: %f; mean: %f; median: %f' %(np.nanmin(tau2_dat), np.nanmax(tau2_dat), np.nanmean(tau2_dat), np.nanmedian(tau2_dat)))

ft.writeto(fitsfile1+'_opdep.fits', tau1_dat, hdr1, output_verify='fix+warn', overwrite=True)
ft.writeto(fitsfile2+'_opdep.fits', tau2_dat, hdr2, output_verify='fix+warn', overwrite=True)
