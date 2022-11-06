# -*- coding: utf-8 -*-
#
# Copyright 2022 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Compute intersection of DR3 query with the 200M star sample
"""

from astropy.table import join, setdiff, Column
import numpy as np


def step4(r_dr3, r_200m, verbose):
    """Compute intersection of DR3 query with the 200M star sample

    Parameters
    ----------
    r_dr3 : astropy Table
        Table containing the initial EDR3 query result.
    r_200m : astropy Table
        Table containing the objects belonging to the 200M star sample.
    verbose : bool
        If True, display additional information.

    Returns
    -------
    r_dr3_200m : astropy Table
        Objects in both the DR3 query and the 200M sample.
    r_dr3_no200m : astropy Table
        Objects in the DR3 query no present in the 200M sample.

    """
    print('<STEP4> Cross-matching DR3 with 200M sample')

    r_200m.rename_column('sourceid', 'source_id')

    # objects in DR3 query also available in the 200M sample
    r_dr3_200m = join(r_dr3, r_200m, keys='source_id', join_type='inner')

    # all objects in DR3 query with data from the 200M sample when available,
    # or with missing entries when no present in the 200M sample
    r_dr3_all = join(r_dr3, r_200m, keys='source_id', join_type='left')

    # objects in DR3 query no present in the 200M sample
    r_dr3_no200m = setdiff(r_dr3_all, r_dr3_200m, keys=['source_id'])

    print(f'        --> Number of objects in the 200M subsample............: {len(r_200m)}')
    print(f'        --> Number of objects in DR3 query.....................: {len(r_dr3_all)}')
    print(f'        --> Number of DR3 objects within the 200M sample.......: {len(r_dr3_200m)}')
    print(f'        --> Number of DR3 objects no present in the 200M sample: {len(r_dr3_no200m)}')

    # sort tables by RA
    r_dr3_all.sort('ra')
    r_dr3_200m.sort('ra')
    r_dr3_no200m.sort('ra')

    # include first column with sequential number (useful for the plot)
    r_dr3_200m.add_column(
        Column(np.arange(1, len(r_dr3_200m)+1)),
        name='number',
        index=0
    )
    r_dr3_no200m.add_column(
        Column(np.arange(1, len(r_dr3_no200m)+1)),
        name='number',
        index=0
    )

    if verbose:
        print('\n* Displaying r_dr3_all:')
        r_dr3_all.pprint()
        print('\n* Displaying r_dr3_200m:')
        r_dr3_200m.pprint()
        print('\n* Displaying r_dr3_no200m:')
        r_dr3_no200m.pprint()

    return r_dr3_200m, r_dr3_no200m
