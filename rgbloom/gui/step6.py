# -*- coding: utf-8 -*-
#
# Copyright 2022 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Generate output PDF plot
"""

from astropy import units as u
from astropy.table import vstack
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style
import numpy as np

from ..choices_mag_plot import CHOICES_MAG_PLOT
from .style import mpl_style
OUTTYPES_COLOR = {'200m': 'red', 'no200m': 'black', 'var': 'blue'}

matplotlib.use('pdf')


def step6(r_dr3_200m, r_dr3_no200m, ra_center, dec_center, search_radius, brightlimit,
          symbsize, max_symbsize, min_symbsize, mag_power,
          display_mag, num_fontsize, nonumbers, nocolor,
          basename, version, verbose):
    """Perform EDR3 query

    Parameters
    ----------
    r_dr3_200m : astropy Table
        Objects in both the DR3 query and the 200M sample.
    r_dr3_no200m : astropy Table
        Objects in the DR3 query no present in the 200M sample.
    ra_center : float
        Right ascension (decimal degree) corresponding to the center
        of the field of view.
    dec_center : float
        Declination (decimal degree) corresponding to the center
        of the field of view.
    search_radius : float
        Radius (decimal degrees) of the field of view.
    brightlimit : float
        Stars brighter than this Gaia G limit are displayed with star
        symbols.
    symbsize : float
        Global multiplying factor for symbol size.
    max_symbsize : float
        Maximum symbol size in chart.
    min_symbsize : float
        Minimum symbol size in chart.
    mag_power: float
        Power to scale symbol sizes in chart. The relative magnitude
        difference (rescaled between 0 and 1) is raised to this power
        to set the symbol size.
    display_mag : bool
        If True display selected magnitude instead of object number.
    num_fontsize : int
        Font size for numbers in chart.
    nonumbers : bool
        If True, do not display star numbers in PDF chart.
    nocolor : bool
        If True, do not use colors in PDF chart.
    basename : str
        Base name for output files.
    version : str
        Version number.
    verbose : bool
        If True, display additional information.

    """
    print('<STEP6> Generating PDF plot')

    # color map
    cmap = plt.cm.get_cmap('jet')
    cmap_vmin = -0.5
    cmap_vmax = 2.0

    # define WCS
    naxis1 = 1024
    naxis2 = naxis1
    pixscale = 2 * search_radius / naxis1

    wcs_image = WCS(naxis=2)
    wcs_image.wcs.crpix = [naxis1 / 2, naxis2 / 2]
    wcs_image.wcs.crval = [ra_center, dec_center]
    wcs_image.wcs.cunit = ["deg", "deg"]
    wcs_image.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs_image.wcs.cdelt = [-pixscale, pixscale]
    wcs_image.array_shape = [naxis1, naxis2]
    if verbose:
        print(wcs_image)

    matplotlib.style.use(mpl_style)
    fig = plt.figure(figsize=(13, 10))
    ax = plt.subplot(projection=wcs_image)

    # generate plot
    if len(r_dr3_200m) > 0 and len(r_dr3_no200m) > 0:
        r_table = vstack([r_dr3_200m, r_dr3_no200m])
    elif len(r_dr3_200m) > 0:
        r_table = r_dr3_200m
    elif len(r_dr3_no200m) > 0:
        r_table = r_dr3_no200m
    else:
        r_table = None

    if r_table is not None:
        r_table.sort('phot_g_mean_mag')
        if verbose:
            r_table.pprint(max_width=1000)

        # define symbol size
        min_mag = np.ma.min(r_table['phot_g_mean_mag'].value)
        max_mag = np.ma.max(r_table['phot_g_mean_mag'].value)
        norm_delta_mag = (np.array(r_table['phot_g_mean_mag'])-min_mag)/(max_mag - min_mag)
        symbol_size = max_symbsize - (norm_delta_mag**(1/mag_power)) * (max_symbsize - min_symbsize)
        symbol_size *= symbsize

        # (X, Y) coordinates
        ra_array = np.array(r_table['ra'])
        dec_array = np.array(r_table['dec'])
        c = SkyCoord(ra=ra_array * u.degree, dec=dec_array * u.degree, frame='icrs')
        x_pix, y_pix = wcs_image.world_to_pixel(c)

        iok_bool = r_table['phot_g_mean_mag'] > brightlimit
        iok = np.arange(len(iok_bool))[iok_bool]
        iokstar = np.arange(len(iok_bool))[~iok_bool]
        if nocolor:
            for i in iokstar:
                ax.scatter(x_pix[i], y_pix[i], marker='*', color='grey',
                           edgecolors='white', linewidth=0.2, s=symbol_size[i], zorder=i+1)
            for i in iok:
                ax.scatter(x_pix[i], y_pix[i], marker='.', color='grey',
                           edgecolors='white', linewidth=0.2, s=symbol_size[i], zorder=i+1)
        else:
            for i in iokstar:
                ax.scatter(x_pix[i], y_pix[i], marker='*',
                           edgecolors='black', linewidth=0.2, s=symbol_size[i], zorder=i+1,
                           cmap=cmap, c=r_table[i]['bp_rp'], vmin=cmap_vmin, vmax=cmap_vmax)
            for i in iok:
                ax.scatter(x_pix[i], y_pix[i], marker='.',
                           edgecolors='black', linewidth=0.2, s=symbol_size[i], zorder=i+1,
                           cmap=cmap, c=r_table[i]['bp_rp'], vmin=cmap_vmin, vmax=cmap_vmax)

        # display numbers if requested
        if not nonumbers:
            for irow in range(len(r_table)):
                if CHOICES_MAG_PLOT[display_mag] is not None:
                    dumvalue = r_table[irow][CHOICES_MAG_PLOT[display_mag]]
                    if isinstance(dumvalue, np.ma.core.MaskedConstant):
                        text = 'NA'  # Not Available
                    else:
                        text = f'{dumvalue:.2f}'
                else:
                    text = r_table[irow]['number']
                if r_table[irow]['qlflag'] == 1:
                    bbox = dict(facecolor='none', edgecolor='gray', boxstyle='round, pad=0.2', lw=1, alpha=0.3)
                else:
                    bbox = None
                if isinstance(r_table[irow]['qlflag'], np.int64):
                    textcolor = OUTTYPES_COLOR['200m']
                else:
                    textcolor = OUTTYPES_COLOR['no200m']
                ax.text(x_pix[irow], y_pix[irow], text, bbox=bbox,
                        color=textcolor, fontsize=num_fontsize,
                        horizontalalignment='left', verticalalignment='bottom',
                        zorder=len(r_table) + 10)

    if len(r_dr3_no200m) > 0:
        # (X, Y) coordinates
        ra_array = np.array(r_dr3_no200m['ra'])
        dec_array = np.array(r_dr3_no200m['dec'])
        c = SkyCoord(ra=ra_array * u.degree, dec=dec_array * u.degree, frame='icrs')
        x_pix, y_pix = wcs_image.world_to_pixel(c)

        # stars outside the -0.5 < G_BP - G_RP < 2.0 colour cut
        mask_colour = np.logical_or((r_dr3_no200m['bp_rp'] <= -0.5), (r_dr3_no200m['bp_rp'] >= 2.0))
        if np.any(mask_colour):
            iok_bool = np.argwhere(mask_colour)
            ax.scatter(x_pix[iok_bool], y_pix[iok_bool], s=240, marker='D',
                       facecolors='none', edgecolors='grey', linewidth=0.5, zorder=0)

        # variable stars
        mask_variable = r_dr3_no200m['phot_variable_flag'] == 'VARIABLE'
        if np.any(mask_variable):
            iok_bool = np.argwhere(mask_variable)
            ax.scatter(x_pix[iok_bool], y_pix[iok_bool], s=240, marker='s',
                       facecolors='none', edgecolors=OUTTYPES_COLOR['var'], linewidth=0.5, zorder=0)

    # legend
    ax.scatter(0.03, 0.97, s=240, marker='s', facecolors='white',
               edgecolors=OUTTYPES_COLOR['var'], linewidth=0.5,
               transform=ax.transAxes)
    ax.text(0.06, 0.97, 'variable in Gaia DR3', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    ax.scatter(0.03, 0.93, s=240, marker='D', facecolors='white', edgecolors='grey', linewidth=0.5,
               transform=ax.transAxes)
    ax.text(0.06, 0.93, 'outside colour range', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    if CHOICES_MAG_PLOT[display_mag] is None:
        ax.text(0.03, 0.89, 'n', color='red', fontsize=12,
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.06, 0.89, '# in *_200m.csv', color='gray', fontsize=12,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.03, 0.85, 'n', color='black', fontsize=12,
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.06, 0.85, '# in *_no200m.csv', color='gray', fontsize=12,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
    else:
        ax.text(0.03, 0.89, 'mag', color='red', fontsize=12,
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.06, 0.89, f'{display_mag} in *_200m.csv', color='gray', fontsize=12,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.03, 0.85, 'mag', color='black', fontsize=12,
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.06, 0.85, f'{display_mag} in *_no200m.csv', color='gray', fontsize=12,
                horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    # plot labels
    ax.set_xlabel('ra')
    ax.set_ylabel('dec')

    ax.set_aspect('equal')

    if not nocolor:
        cbaxes = fig.add_axes([0.683, 0.81, 0.15, 0.02])
        norm = mpl.colors.Normalize(vmin=cmap_vmin, vmax=cmap_vmax)
        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbaxes,
                            orientation='horizontal', format='%1.0f')
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(label=r'$G_{\rm BP}-G_{\rm RP}$ (mag)', size=12, backgroundcolor='white')

    ax.text(0.98, 0.96, f'Field radius: {search_radius:.4f} degree', fontsize=12, backgroundcolor='white',
            horizontalalignment='right', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.02, 0.06, r'$\alpha_{\rm center}$:', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    ax.text(0.25, 0.06, f'{ra_center:.4f} degree', fontsize=12, backgroundcolor='white',
            horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
    ax.text(0.02, 0.02, r'$\delta_{\rm center}$:', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    ax.text(0.25, 0.02, f'{dec_center:+.4f} degree', fontsize=12, backgroundcolor='white',
            horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
    ax.text(0.98, 0.02, f'rgbloom, version {version}', fontsize=12, backgroundcolor='white',
            horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)

    f = np.pi / 180
    xp = naxis1 / 2 + search_radius/pixscale * np.cos(np.arange(361)*f)
    yp = naxis2 / 2 + search_radius/pixscale * np.sin(np.arange(361)*f)
    ax.plot(xp, yp, '-', color='orange', linewidth=0.5, alpha=0.5)

    ax.set_xlim([-naxis1*0.12, naxis1*1.12])
    ax.set_ylim([-naxis2*0.05, naxis2*1.05])

    ax.set_axisbelow(True)
    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='black', ls='dotted')

    plt.savefig(f'{basename}.pdf')
    plt.close(fig)
