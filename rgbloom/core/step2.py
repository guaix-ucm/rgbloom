# -*- coding: utf-8 -*-
#
# Copyright 2022 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Estimate RGB magnitudes using C21 polynomials
"""

from astropy import units as u
from astropy.table import Column
import numpy as np
from numpy.polynomial import Polynomial
import sys


def step2(r_dr3, verbose):
    """Estimate RGB magnitudes using C21 polynomials

    Parameters
    ----------
    r_dr3 : astropy Table
        Table containing the initial DR3 query.
    verbose : bool
        If True, display additional information.

    Returns
    -------
    r_dr3 : astropy Table
        Updated DR3 table including the predicted RGB magnitudes.

    """
    sys.stdout.write('<STEP2> Estimating RGB magnitudes in DR3 query using C21 polynomials')
    sys.stdout.flush()
    # predict RGB magnitudes
    coef_b = np.array([-0.13748689, 0.44265552, 0.37878846, -0.14923841, 0.09172474, -0.02594726])
    coef_g = np.array([-0.02330159, 0.12884074, 0.22149167, -0.1455048, 0.10635149, -0.0236399])
    coef_r = np.array([0.10979647, -0.14579334, 0.10747392, -0.1063592, 0.08494556, -0.01368962])

    poly_b = Polynomial(coef_b)
    poly_g = Polynomial(coef_g)
    poly_r = Polynomial(coef_r)

    r_dr3.add_column(
        Column(np.round(r_dr3['phot_g_mean_mag'] + poly_b(r_dr3['bp_rp']), 3),
               name='RGB_B_C21', unit=u.mag, format='.3f')
    )
    r_dr3.add_column(
        Column(np.round(r_dr3['phot_g_mean_mag'] + poly_g(r_dr3['bp_rp']), 3),
               name='RGB_G_C21', unit=u.mag, format='.3f')
    )
    r_dr3.add_column(
        Column(np.round(r_dr3['phot_g_mean_mag'] + poly_r(r_dr3['bp_rp']), 3),
               name='RGB_R_C21', unit=u.mag, format='.3f')
    )
    print(' OK!')
    if verbose:
        r_dr3.pprint(max_width=1000)

    return r_dr3
