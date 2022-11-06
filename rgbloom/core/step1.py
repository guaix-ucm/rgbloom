# -*- coding: utf-8 -*-
#
# Copyright 2022 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Perform DR3 query
"""

import sys

from astropy import units as u
from astropy.table import Column
from astroquery.gaia import Gaia
import numpy as np


def step1(ra_center, dec_center, search_radius, g_limit, verbose):
    """Perform DR3 query

    Parameters
    ----------
    ra_center : float
        Right ascension (decimal degree) corresponding to the center
        of the field of view.
    dec_center : float
        Declination (decimal degree) corresponding to the center
        of the field of view.
    search_radius : float
        Radius (decimal degrees) of the field of view.
    g_limit : float
        Limiting Gaia G magnitude.
    verbose : bool
        If True, display additional information.

    Returns
    -------
    r_dr3 : astropy Table
        Table containing the query result.
    nstars : int
        Number of stars in 'r_dr3'.

    """
    query = f"""
    SELECT source_id, ra, dec,
    phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
    phot_variable_flag

    FROM gaiadr3.gaia_source_lite
    WHERE 1=CONTAINS(
      POINT('ICRS', {ra_center}, {dec_center}), 
      CIRCLE('ICRS',ra, dec, {search_radius}))
    AND phot_g_mean_mag IS NOT NULL 
    AND phot_bp_mean_mag IS NOT NULL 
    AND phot_rp_mean_mag IS NOT NULL
    AND phot_g_mean_mag < {g_limit}

    ORDER BY ra
    """
    sys.stdout.write('<STEP1> Starting cone search in Gaia DR3... (please wait)\n  ')
    sys.stdout.flush()
    job = Gaia.launch_job_async(query)
    r_dr3 = job.get_results()
    nstars = len(r_dr3)
    if nstars == 0:
        raise SystemExit('ERROR: no objects found. Change search parameters!')
    # compute G_BP - G_RP colour
    r_dr3.add_column(
        Column(r_dr3['phot_bp_mean_mag'] - r_dr3['phot_rp_mean_mag'],
               name='bp_rp', unit=u.mag)
    )
    # colour cut in BP-RP
    mask_colour = np.logical_or((r_dr3['bp_rp'] <= -0.5), (r_dr3['bp_rp'] >= 2.0))
    r_dr3_colorcut = r_dr3[mask_colour]
    print(f'        --> {nstars} objects found')
    # check for variable objects
    mask_variable = r_dr3['phot_variable_flag'] == 'VARIABLE'
    r_dr3_variable = r_dr3[mask_variable]
    nstars_variable = len(r_dr3_variable)
    print(f'        --> {nstars_variable} objects classified as VARIABLE')
    if verbose:
        r_dr3.pprint(max_width=1000)

    return r_dr3, nstars
