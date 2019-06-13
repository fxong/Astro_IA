"""
Purpose:
    plot the radial density profiles of filaments
Usage:
Caution:
History:
    Completed. Version 0.9.0. Feb. 25, 19.
Copyright:
    written by fxong@CfA
"""

import os
pwd=os.getcwd()+'/'
os.chdir(pwd)
import numpy as np
import matplotlib as mpl
mpl.rcParams['xtick.direction']='in'
mpl.rcParams['ytick.direction']='in'
import matplotlib.pyplot as plt

from astropy.io import ascii
from astropy.io import fits as ft
from astropy.wcs import WCS as wcs

tablefile='FilProfile_profile.data' # table of the profile of filament
int_range=[] # integrated range of line mass, unit is pc, default is []
col_thres=0 # column density threshold of line mass, unit is cm^-2
figsize=None # size of the figure, e.g., [6.4, 4.8]/[8.0, 6.0], default is None

""" inputs """
fil_path=ascii.read(tablefile, format='fixed_width', names=['x', 'y', 'stdev', 'err'])
x=np.array([i for i in fil_path['x']])
y=np.array([i for i in fil_path['y']])
stdev=np.array([i for i in fil_path['stdev']])
err=np.array([i for i in fil_path['err']])

if int_range==[]:
    maxvalue=np.nanmax(y[(x > -0.1) & (x < 0.1)])
    negin=np.where(y==maxvalue)[0]
    while y[negin]>col_thres: negin=negin-1
    posin=np.where(y==maxvalue)[0]
    while y[posin]>col_thres: posin=posin+1
    int_range=[x[negin], x[posin]]
# elif:
#     negin=np.nanargmin(abs(x-int_range[0]))
#     posin=np.nanargmin(abs(x-int_range[1]))

""" line mass """
int_x=x[(x >= int_range[0]) & (x <= int_range[1])]
int_y=y[(x >= int_range[0]) & (x <= int_range[1])]
int_err=err[(x >= int_range[0]) & (x <= int_range[1])]

col_y=(int_y*2.83*(1.6726e-27)/(1.9885e30))/((1/3.0856776e18)**2)
col_err=(int_err*2.83*(1.6726e-27)/(1.9885e30))/((1/3.0856776e18)**2)

lin_ms=0
lin_err=0
for i in range(int_x.shape[0]-1):
    lin_ms=lin_ms+np.abs(int_x[i+1]-int_x[i])*(col_y[i+1]+col_y[i])/2
    lin_err=lin_err+np.power(np.abs(int_x[i+1]-int_x[i]),2)*(np.power(col_err[i+1],2)+np.power(col_err[i],2))/4
lin_err=np.sqrt(lin_err)

print('\nThe integrated range of line mass is from %f pc to %f pc' %(int_range[0], int_range[1]))
print('The line mass of the filament is %f +/- %f Msun/pc\n' %(lin_ms, lin_err))

""" import figure and subplot """
fig=plt.figure(figsize=figsize)
ax=plt.subplot(1, 1, 1)
ax.set_aspect('auto')

xdat=x
ydat=y/1e23
stdevy=stdev/1e23

""" axis labels and tick labels """
xlabel='Radius (pc)'
ylabel='$N_{\\rm H_2}$ (10$^{23}$ cm$^{-2}$)'
ax.set_xlabel(xlabel, labelpad=None, size='medium', weight='normal', visible=True)
ax.set_ylabel(ylabel, labelpad=None, size='medium', weight='normal', visible=True)

ax.tick_params(axis='x', labelsize='medium', labelbottom=True)
ax.tick_params(axis='y', labelsize='medium', labelleft=True)

""" scaling and ticks """
ax.set_xscale(value='linear')
ax.set_yscale(value='linear')
ax.set_xlim(-0.8, 0.8)
ax.set_ylim(-0.02, 4.3)

ax.tick_params(axis='x', which='both', top=True)
ax.tick_params(axis='y', which='both', right=True)
# ax.tick_params(axis='x', which='major', length=, width=)
# ax.tick_params(axis='y', which='major', length=, width=)
ax.minorticks_on()

""" scatter and plot """
ax.errorbar(xdat, ydat, fmt='none', yerr=stdevy, ecolor='gold', elinewidth=1.2, capsize=2.0, capthick=1.2)

ax.plot(xdat, ydat, color='black', linestyle='solid', linewidth=1.5, label='Observed profile')

ax.plot([int_range[0], int_range[0]], ax.get_ylim(), color='black', linestyle='dashed', linewidth=1.0, label='')
ax.plot([int_range[1], int_range[1]], ax.get_ylim(), color='black', linestyle='dashed', linewidth=1.0, label='')

""" model plotting """
Npl=
p=
Rf=
Plx=np.linspace(np.nanmin(xdat), np.nanmax(xdat), xdat.shape[0]*3)
Ply=Npl/np.power(1+(Plx/Rf)**2, (p-1)/2)

ax.plot(Plx, Ply, color='red', linestyle='dashed', linewidth=1.5, label='Plummer fit')

Ngs=
mu=
sigma=
Gsx=np.linspace(np.nanmin(xdat), np.nanmax(xdat), xdat.shape[0]*3)
Gsy=Ngs*np.exp((-(Gsx-mu)**2)/(2*(sigma**2)))

ax.plot(Gsx, Gsy, color='dodgerblue', linestyle='dashed', linewidth=1.5, label='Gaussian fit')

ax.legend(fontsize='medium', facecolor=None, edgecolor=None, framealpha=0.0, loc='upper right', borderaxespad=0.2)

""" export figure file """
plt.show()

os.system('rm -f FilPlot.pdf')
fig.savefig('FilPlot.pdf', format='pdf')
