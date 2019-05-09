"""
Purpose:
    plot the input 3D fitsfile and export a figure file
Usage:
Caution:
    fitsfile must be 3D
    most arguments need to be set in the command
History:
    Completed. Version 0.8.0. Nov. 05, 18.
    Updated. Version 0.9.0. Nov. 30, 18.
Copyright:
    written by fxong@CfA
"""

import os
pwd=os.getcwd()+'/'
os.chdir(pwd)

import numpy as np
import astropy.units as u
from astropy.io import ascii
from astropy.io import fits as ft
from astropy.wcs import WCS as wcs
from astropy.wcs.utils import proj_plane_pixel_scales as pixsc

import matplotlib as mpl
mpl.rcParams['xtick.direction']='in'
mpl.rcParams['ytick.direction']='in'
import matplotlib.pyplot as plt

from astropy.visualization import ImageNormalize
from astropy.visualization import ManualInterval
from astropy.visualization import LinearStretch, LogStretch, PowerStretch, SqrtStretch, SquaredStretch

from matplotlib.colorbar import make_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse as AE
# from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDirectionArrows as ADA

from scipy.ndimage import gaussian_filter

fitsfile=''
# fitsfile without suffix
figsize=[8.0, 6.0]
# size of the figure, default is [8.0, 6.0]
scale=None
# scale of the datacube when plotting, e.g., [min, max], default is None
colormap=None
# colormap when plotting, e.g., 'rainbow', default is None
row_col=[3, 3]
# rows and columns of the panels in the figure, default is [3, 3]
outfig=fitsfile
# output eps file without suffix, default is fitsfile.pdf

dat, hdr=ft.getdata(fitsfile+'.fits', header=True)
hdr['NAXIS']=2
del hdr['NAXIS3'], hdr['PC3_1'], hdr['PC3_2'], hdr['PC1_3'], hdr['PC2_3'], hdr['PC3_3']
del hdr['CTYPE3'], hdr['CRVAL3'], hdr['CDELT3'], hdr['CRPIX3'], hdr['CUNIT3']

distance=1*1e3*u.pc
# distance to the source, unit is 1e3*u.pc
pixsz=np.mean(pixsc(wcs(hdr)))*u.deg
# pixel size of the datacube, unit is u.deg

""" import figure and subplot """
fig=plt.figure(figsize=figsize)
for pltid in range(dat.shape[0]):
    ax=plt.subplot(row_col[0], row_col[1], pltid+1, projection=wcs(hdr))
    plt.subplots_adjust(wspace=0, hspace=0)
    xaxis, yaxis, xyaxis=ax.coords[0], ax.coords[1], ax.coords

    """ axis labels and tick labels """
    show=True if pltid==(row_col[0]-1)*row_col[1] else False
    xax_vis, yax_vis=show, show
    # whether to display the axis labels or not, defaut is True
    xtk_vis, ytk_vis=show, show
    # whether to display the tick labels or not, defaut is True
    xpad, ypad=1.0, 1.0
    # distance between axis labels and tick labels, defaut is 1.0
    xsize, ysize='medium', 'medium'
    # size of the labels, e.g., 'x-small'/'small'/'medium'/'large'/'x-large', defaut is 'medium'
    xweit, yweit='normal', 'normal'
    # weight of the labels, e.g., 'ultralight'/'light'/'normal'/'regular'/'book', defaut is 'normal'
    ytk_rot=0
    # rotation of y axis tick labels, unit is arcdeg, defaut is 0

    xlabel='Right Ascension (J2000)' # Right Ascension (J2000)
    ylabel='Declination (J2000)' # Declination (J2000)
    xaxis.set_axislabel(xlabel, minpad=xpad, size=xsize, weight=xweit, visible=xax_vis)
    yaxis.set_axislabel(ylabel, minpad=ypad, size=ysize, weight=yweit, visible=yax_vis)

    xaxis.set_major_formatter('hh:mm:ss.s')
    yaxis.set_major_formatter('dd:mm:ss.s')
    # xaxis.set_separator(('$^{\\rm h}$', '$^{\\rm m}$', '$^{\\rm s}$'))
    xaxis.set_ticklabel(size=xsize, weight=xweit, visible=xtk_vis)
    yaxis.set_ticklabel(size=ysize, weight=yweit, visible=ytk_vis, rotation=ytk_rot)

    """ major ticks and minor ticks """
    xma_num, yma_num=None, None
    # number of major ticks, default is None
    xcolor, ycolor='black', 'black'
    # color of major and minor ticks, default is 'black'
    xlenth, ylenth=None, None
    # length of major and minor ticks, default is None
    xwidth, ywidth=None, None
    # width of major and minor ticks, default is None
    xmi_num, ymi_num=None, None
    # number of minor ticks, default is None

    xaxis.set_ticks(number=xma_num, color=xcolor, size=xlenth, width=xwidth, exclude_overlapping=True)
    yaxis.set_ticks(number=yma_num, color=ycolor, size=ylenth, width=ywidth, exclude_overlapping=True)

    xaxis.display_minor_ticks(True)
    yaxis.display_minor_ticks(True)
    # xaxis.set_minor_frequency(xmi_num)
    # yaxis.set_minor_frequency(ymi_num)

    """ scaling and imaging """
    interp='none'
    # interpolation of the datacube, e.g., 'none'/'nearest'/'bilinear'/'bicubic', defaut is 'none'

    scaling='linear'
    # scaling of the datacube, e.g., 'linear'/'log'/'power'/'sqrt'/'squared', defaut is 'linear'
    index=1e3
    # index of log and power scaling, default is 1e3
    if scaling=='linear': stretch=LinearStretch()
    elif scaling=='log': stretch=LogStretch(index)
    elif scaling=='power': stretch=PowerStretch(index)
    elif scaling=='sqrt': stretch=SqrtStretch()
    elif scaling=='squared': stretch=SquaredStretch()
    if scale is None: interval=ManualInterval(vmin=np.nanmin(dat), vmax=np.nanmax(dat))
    else: interval=ManualInterval(vmin=scale[0], vmax=scale[1])
    norm=ImageNormalize(dat, interval=interval, stretch=stretch)

    # ax.set_xlim(-0.5, dat.shape[1]-0.5)
    # ax.set_ylim(-0.5, dat.shape[0]-0.5)
    if colormap is not None: plt.set_cmap(colormap)
    cmap=mpl.cm.get_cmap()
    cmap.set_bad(color='white')
    im=ax.imshow(dat, norm=norm, interpolation=interp, cmap=cmap, origin='lower')

    """
    The following arguments need to be set in the command.
    """

    """ colorbar and beam """
    cbpos='right' # position of colorbar, 'right'/'top'
    cbaxpos=0.9 # modifying position of colorbar, default is 0.9
    dlenth=0 # delta length of colorbar, default is 0
    cbwidth=0.05 # width of colorbar, default is 0.05
    cblab='Jy beam$^{-1}$ km s$^{-1}$' # label of colorbar
    ticks=None # tick labels of colorbar, default is None
    # cax.minorticks_on() # whether to show the minor ticks or not

    if pltid==row_col[1]-1:
        if cbpos=='right':
            orientation, cbx_vis, cby_vis='vertical', False, True
            cax, kw=make_axes(ax, location='right', pad=0.0)
            cbox=cax.get_position(True)
            cbox.x0, cbox.x1=cbaxpos, cbaxpos+cbox.x1-cbox.x0
            cbox.y0, cbox.y1=cbox.y0+(dlenth/2.0), cbox.y1-(dlenth/2.0)
            aspcet=1/cbwidth
        elif cbpos=='top':
            orientation, cbx_vis, cby_vis='horizontal', True, False
            cax, kw=make_axes(ax, location='top', pad=0.0)
            cbox=cax.get_position(True)
            cbox.x0, cbox.x1=cbox.x0+(dlenth/2.0), cbox.x1-(dlenth/2.0)
            cbox.y0, cbox.y1=cbaxpos, cbaxpos+cbox.y1-cbox.y0
            aspcet=cbwidth
        cax.set_position(cbox)
        cax.set_aspect(aspcet)
        plt.colorbar(im, cax=cax, orientation=orientation, ticks=ticks)

        cax.xaxis.set_label_position('top')
        cax.set_xlabel(cblab, size='medium', weight='normal', labelpad=None, visible=cbx_vis)
        cax.yaxis.set_label_position('right')
        cax.set_ylabel(cblab, size='medium', weight='normal', labelpad=None, visible=cby_vis)

        # if cbpos=='right': ticks=cax.get_yticklabels()
        # elif cbpos=='top': ticks=cax.get_xticklabels()
        # ticks=[float(i.get_text()) for i in ticks]
        # print('The tick labels of colorbar are', ticks)
        cax.tick_params(axis='x', which='both', color='black', top=True, bottom=False, labelsize='medium', labeltop=True, labelbottom=False, labelrotation=0)
        cax.tick_params(axis='y', which='both', color='black', right=True, left=False, labelsize='medium', labelright=True, labelleft=False, labelrotation=0)

    major_fwhm=hdr['bmaj']*u.deg # hdr['bmaj']*u.deg
    minor_fwhm=hdr['bmin']*u.deg # hdr['bmin']*u.deg
    angle=hdr['bpa']*u.deg # hdr['bpa']*u.deg

    beam=AE(ax.transData, height=major_fwhm/pixsz, width=minor_fwhm/pixsz, angle=angle, loc='lower left', pad=0.2, borderpad=0.2, frameon=True)
    beam.ellipse.set(color='black', hatch='/////', linewidth=1.0, fill=True)
    beam.patch.set(edgecolor='black', linewidth=1.0, fill=False)
    ax.add_artist(beam)

    """ scalebar and contour """
    barlen=0*u.pc # u.pc/u.deg/u.arcmin/u.arcsec
    # bardeg=((barlen/distance)/np.pi)*(180*u.deg)
    # bardeg=barlen
    # bardeg=(barlen/(60*u.arcmin))*u.deg
    # bardeg=(barlen/(3600*u.arcsec))*u.deg
    barlab=str(barlen.value)+' pc' # ' pc'/'$^{\\circ}$'/"'"/'"'

    # ax.hlines(y, x, x+(bardeg/pixsz), colors='black', linestyles='solid', linewidths=1.0, transform=ax.get_transform('pixel'))
    # ax.vlines(x, y, y+(bardeg/pixsz), colors='black', linestyles='solid', linewidths=1.0, transform=ax.get_transform('pixel'))
    # ax.text(x, y, barlab, color='black', size='medium', weight='normal', rotation=0, transform=ax.get_transform('pixel'))

    levels=[]
    contfile=''
    smooth=0.0 # smoothing level of the contours (under development!!), set to be <1.0, defaut is 0.0

    # cdat, chdr=ft.getdata(contfile+'.fits', header=True)
    # cdat=gaussian_filter(cdat, sigma=smooth)
    # ax.contour(cdat, colors=['black'], linewidths=[1.0], linestyles=['solid'])
    # ax.contour(cdat, levels=levels, colors=['black'], linewidths=[1.0], linestyles=['solid'])
    # ax.contour(cdat, levels=levels, colors=['black'], linewidths=[1.0], linestyles=['solid'], transform=ax.get_transform(wcs(chdr)))

    """ markers """
    # xyaxis.grid(color='black', linestyle='dashed', linewidth=1.0)

    #tabxy=ascii.read('xy.table', format='fixed_width', names=['x', 'y'])

    # ax.scatter(x, y, s=36, marker='o', facecolor='black', edgecolor='none', linewidths=1.0, label='', transform=ax.get_transform('icrs'))

    # ax.plot(x, y, color='black', linestyle='solid', linewidth=1.0, label='', transform=ax.get_transform('icrs'))

    # ax.text(x, y, 'text', color='black', size='medium', weight='normal', rotation=0, transform=ax.get_transform('pixel'))

""" export figure file """
# fig.savefig('temp.ps',format='ps')
# os.system('ps2ps2 temp.ps temp.ps2 \n ps2eps temp.ps2')
# os.system('rm -f '+outfig+'.eps')
# os.rename('temp.ps2.eps', outfig+'.eps')
# os.system('rm -f temp.ps temp.ps2')
# plt.show()

os.system('rm -f '+outfig+'.pdf')
fig.savefig(outfig+'.pdf', format='pdf')
