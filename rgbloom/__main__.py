# -*- coding: utf-8 -*-
#
# Copyright 2022-2023 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
RGB predictions of Gaia DR3 stars with XP spectra

This code is hosted at https://github.com/guaix-ucm/rgbloom
Authors: Nicolás Cardiel <cardiel@ucm.es>
         Sergio Pascual <sergiopr@fis.ucm.es>
         Rafael González <rafael08@ucm.es>

Usage example:
$ rgbloom 56.66 24.10 1.0 12
"""

import argparse
import pandas as pd
import pooch
import sys

from .choices_mag_plot import CHOICES_MAG_PLOT
from .core.step1 import step1
from .core.step2 import step2
from .core.step3 import step3
from .core.step4 import step4
from .core.step5 import step5
from .gui.step6 import step6
from .version import version

MAX_SEARCH_RADIUS = 30  # degrees
REFERENCE_HEALPIX8 = "reference_healpix8.csv"


def right_ascension(ra_str):
    ra = float(ra_str)
    if not (0 <= ra <= 360):
        print('Right ascension must be 0 <= ra <= 360 degree')
        msg = f'Right ascension {ra} degree out of range'
        print(msg)
        raise ValueError(msg)
    return ra


def declination(dec_str):
    dec = float(dec_str)
    if not (-90 <= dec <= 90):
        print('Declination must be -90 <= dec <= 360 degree')
        msg = f'Declination {dec} degree out of range'
        print(msg)
        raise ValueError(msg)
    return dec


def search_radius(r_str):
    r = float(r_str)
    if not (0 < r < MAX_SEARCH_RADIUS):
        print(f'Search radius must be 0 < r <= {MAX_SEARCH_RADIUS} degree')
        msg = f'Search radius {r} degree out of range'
        print(msg)
        raise ValueError(msg)
    return r


def exec_rgbloom(args):
    """Callable function that executes all the required steps

    Parameters
    ----------
    args: argparse instance
        Argparse instance containing the command-line parameters

    """
    print(f'\n        Welcome to rgbloom version {version}')
    print(f'        ==============================\n')

    # compute md5 hash from terminal using:
    # linux $ md5sum <filename>
    # macOS $ md5 <filename>
    fauxcsv = pooch.retrieve(
        f"http://nartex.fis.ucm.es/~ncl/rgbphot/gaiaDR3/{REFERENCE_HEALPIX8}",
        known_hash="md5:65dbb1f2c5030a1bc987607fa5783736"
    )
    if args.verbose:
        print(f'- Required file: {fauxcsv}')

    # read the previous file
    try:
        reference_healpix8_df = pd.read_csv(fauxcsv, sep=',', header=0)
    except FileNotFoundError:
        raise SystemExit(f'ERROR: unexpected problem while reading {REFERENCE_HEALPIX8}')

    # ---
    # step 1: perform DR3 query
    r_dr3, nstars = step1(
        args.ra_center,
        args.dec_center,
        args.search_radius,
        args.g_limit,
        args.verbose
    )

    # ---
    # step 2: compute RGB magnitudes in DR3 query
    r_dr3 = step2(
        r_dr3,
        args.verbose
    )

    # ---
    # step 3: retrieve source_id of sources belonging to the 200M sample
    # in the same region of the sky
    r_200m = step3(
        args.ra_center,
        args.dec_center,
        args.search_radius,
        reference_healpix8_df,
        args.verbose
    )

    # ---
    # step 4: compute intersection of DR3 query with the 200M star sample
    r_dr3_200m, r_dr3_no200m = step4(
        r_dr3,
        r_200m,
        args.verbose
    )

    # ---
    # step 5: generate output CSV files
    step5(
        r_dr3_200m,
        r_dr3_no200m,
        args.basename,
        args.verbose,
    )

    # ---
    # step 6: generate PDF chart
    if args.noplot:
        sys.stdout.write('<STEP6> No PDF plot generated (skipped!)\n')
        sys.stdout.flush()
    else:
        step6(
            r_dr3_200m,
            r_dr3_no200m,
            args.ra_center,
            args.dec_center,
            args.search_radius,
            args.brightlimit,
            args.symbsize,
            args.max_symbsize,
            args.min_symbsize,
            args.mag_power,
            args.display_mag,
            args.num_fontsize,
            args.nonumbers,
            args.nocolor,
            args.basename,
            version,
            args.verbose
        )
    print('End of program')


def main():
    """Main function to parse input arguments"""
    parser = argparse.ArgumentParser(description=f"RGB predictions from Gaia DR3 spectrophotometry (version {version})")
    parser.add_argument("ra_center", help="right Ascension (decimal degrees)", type=right_ascension)
    parser.add_argument("dec_center", help="declination (decimal degrees)", type=declination)
    parser.add_argument("search_radius", help="search radius (decimal degrees)", type=search_radius)
    parser.add_argument("g_limit", help="limiting Gaia G magnitude", type=float)
    parser.add_argument("--basename", help="file basename for output files", type=str, default="rgbloom")
    parser.add_argument("--brightlimit",
                        help="objects brighter than this Gaia G limit are displayed with star symbols (default=8.0)",
                        type=float, default=8.0)
    parser.add_argument("--symbsize", help="global multiplying factor for symbol size (default=1.0)",
                        type=float, default=1.0)
    parser.add_argument("--max_symbsize", help="maximum symbol size in chart (default=1000)",
                        type=float, default=1000)
    parser.add_argument("--min_symbsize", help="minimum symbol size in chart (default=10)",
                        type=float, default=10)
    parser.add_argument("--mag_power", help="power to scale symbol sizes in chart (default=3)",
                        type=float, default=3)
    parser.add_argument("--display_mag", help="display selected magnitude instead of object number",
                        type=str, default="None", choices=list(CHOICES_MAG_PLOT.keys()))
    parser.add_argument("--num_fontsize", help="font size for numbers in chart (default=5)",
                        type=int, default=5)
    parser.add_argument("--nonumbers", help="do not display object identification number in PDF chart",
                        action="store_true")
    parser.add_argument("--noplot", help="skip PDF chart generation", action="store_true")
    parser.add_argument("--nocolor", help="do not use colors in PDF chart", action="store_true")
    parser.add_argument("--verbose", help="increase program verbosity", action="store_true")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_usage()
        raise SystemExit()

    exec_rgbloom(args)


if __name__ == "__main__":

    main()
