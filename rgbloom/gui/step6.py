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

import matplotlib.pyplot as plt
import matplotlib.style
from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np

from .style import mpl_style
OUTTYPES_COLOR = {'200m': 'red', 'no200m': 'black', 'var': 'blue'}

matplotlib.use('pdf')


def step6(r_dr3_200m, r_dr3_no200m, ra_center, dec_center, search_radius, symbsize, brightlimit,
          nonumbers, nocolor, basename, version, verbose):
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
    symbsize : float
        Multiplying factor for symbol size.
    brightlimit : float
        Stars brighter than this Gaia G limit are displayed with star
        symbols.
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
    for sample in ['200m', 'no200m']:
        if sample == '200m':
            r_table = r_dr3_200m
        else:
            r_table = r_dr3_no200m

        r_table.sort('phot_g_mean_mag')
        if verbose:
            r_table.pprint(max_width=1000)

        symbol_size = symbsize * (50 / np.array(r_table['phot_g_mean_mag'])) ** 2.5
        ra_array = np.array(r_table['ra'])
        dec_array = np.array(r_table['dec'])

        c = SkyCoord(ra=ra_array * u.degree, dec=dec_array * u.degree, frame='icrs')
        x_pix, y_pix = wcs_image.world_to_pixel(c)

        iok = r_table['phot_g_mean_mag'] < brightlimit
        if nocolor:
            sc = ax.scatter(x_pix[iok], y_pix[iok], marker='*', color='grey',
                            edgecolors='black', linewidth=0.2, s=symbol_size[iok])
            ax.scatter(x_pix[~iok], y_pix[~iok], marker='.', color='grey',
                       edgecolors='black', linewidth=0.2, s=symbol_size[~iok])
        else:
            cmap = plt.cm.get_cmap('jet')
            sc = ax.scatter(x_pix[iok], y_pix[iok], marker='*',
                            edgecolors='black', linewidth=0.2, s=symbol_size[iok],
                            cmap=cmap, c=r_table[iok]['bp_rp'], vmin=-0.5, vmax=2.0)
            ax.scatter(x_pix[~iok], y_pix[~iok], marker='.',
                       edgecolors='black', linewidth=0.2, s=symbol_size[~iok],
                       cmap=cmap, c=r_table[~iok]['bp_rp'], vmin=-0.5, vmax=2.0)

        # display numbers if requested
        if not nonumbers:
            for irow in range(len(r_table)):
                text = r_table[irow]['number']
                if r_table[irow]['qlflag'] == 1:
                    bbox = dict(facecolor='none', edgecolor='gray', boxstyle='round, pad=0.2', lw=1, alpha=0.3)
                else:
                    bbox = None
                ax.text(x_pix[irow], y_pix[irow], text, bbox=bbox,
                        color=OUTTYPES_COLOR[sample], fontsize='5',
                        horizontalalignment='left', verticalalignment='bottom')

        if sample == 'no200m':
            # stars outside the -0.5 < G_BP - G_RP < 2.0 colour cut
            mask_colour = np.logical_or((r_dr3_no200m['bp_rp'] <= -0.5), (r_dr3_no200m['bp_rp'] >= 2.0))
            if np.any(mask_colour):
                iok = np.argwhere(mask_colour)
                ax.scatter(x_pix[iok], y_pix[iok], s=240, marker='D',
                           facecolors='none', edgecolors='grey', linewidth=0.5)

            # variable stars
            mask_variable = r_dr3_no200m['phot_variable_flag'] == 'VARIABLE'
            if np.any(mask_variable):
                iok = np.argwhere(mask_variable)
                ax.scatter(x_pix[iok], y_pix[iok], s=240, marker='s',
                           facecolors='none', edgecolors='blue', linewidth=0.5)

    # legend
    ax.scatter(0.03, 0.96, s=240, marker='s', facecolors='white',
               edgecolors=OUTTYPES_COLOR['var'], linewidth=0.5,
               transform=ax.transAxes)
    ax.text(0.06, 0.96, 'variable in Gaia DR3', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    ax.scatter(0.03, 0.92, s=240, marker='D', facecolors='white', edgecolors='grey', linewidth=0.5,
               transform=ax.transAxes)
    ax.text(0.06, 0.92, 'outside colour range', fontsize=12, backgroundcolor='white',
            horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    ax.set_xlabel('ra')
    ax.set_ylabel('dec')

    ax.set_aspect('equal')

    if not nocolor:
        cbaxes = fig.add_axes([0.683, 0.81, 0.15, 0.02])
        cbar = plt.colorbar(sc, cax=cbaxes, orientation='horizontal', format='%1.0f')
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(label=r'$G_{\rm BP}-G_{\rm RP}$', size=12, backgroundcolor='white')

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
