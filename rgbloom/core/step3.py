# -*- coding: utf-8 -*-
#
# Copyright 2022-2023 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Retrieve sources belonging to the 200M sample
"""

from astropy import units as u
from astropy_healpix import HEALPix
from astropy.table import Table, vstack
import numpy as np
import pandas as pd
import pooch


def step3(ra_center, dec_center, search_radius, reference_healpix8_df, verbose):
    """Retrieve sources belonging to the 200M sample

    The resulting table is constrained to the HEALPIx level-8 tables
    enclosing the same region of the sky.

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
    reference_healpix8_df : pandas dataframe
        Table containing the healpix8_min and healpix8_max values.
    verbose : bool
        If True, display additional information.

    Returns
    -------
    r_200m: astropy Table
        Table containing the objects belonging to the 200M sample in
        the HEALPIx level-8 tables enclosing the same region of the
        sky.
    """
    hpx_level = 8
    hp = HEALPix(nside=2 ** hpx_level, order='nested')
    hp_cone_search = hp.cone_search_lonlat(ra_center * u.deg, dec_center * u.deg, search_radius * u.deg)
    hp_cone_search.sort()

    subset_list = []
    md5_list = []
    for index in reference_healpix8_df.index:
        row = reference_healpix8_df.iloc[index]
        hp_min, hp_max, md5_hash = row['healpix8_min'], row['healpix8_max'], row['md5_hash']
        if np.any(np.logical_and(hp_min <= hp_cone_search, hp_cone_search <= hp_max)):
            bulk_file = f'{hp_min:06d}-{hp_max:06d}'
            subset_list.append(bulk_file)
            md5_list.append(md5_hash)

    print('<STEP3> Retrieving objects from the 200M sample in the enclosing HEALPIx level-8 tables')
    r_200m = None
    subdir = 'RGBsynthetic_NOVARIABLES'
    for subset, md5 in zip(subset_list, md5_list):
        # retrieve file containing the RGB predictions for objects
        # within the 200M sample in a particular region of the sky
        fname = f'{subdir}/sortida_XpContinuousMeanSpectrum_{subset}_RGB_NOVARIABLES_final.csv.gz'
        fauxcsv = pooch.retrieve(
            f"http://nartex.fis.ucm.es/~ncl/rgbphot/gaiaDR3/{fname}",
            known_hash=f"md5:{md5}"
        )
        print(f'        * Required file: {fauxcsv}')
        print(f'          md5:{md5}')
        df = pd.read_csv(fauxcsv, compression='gzip', sep=',', header=0)
        print(f'        --> Number of objects: {len(df)}')
        if r_200m is None:
            r_200m = Table.from_pandas(df)
        else:
            r_200m = vstack([r_200m, Table.from_pandas(df)], join_type='exact')

    print(f'        --> Total number of objects: {len(r_200m)}')
    if verbose:
        r_200m.pprint(max_width=1000)

    return r_200m
