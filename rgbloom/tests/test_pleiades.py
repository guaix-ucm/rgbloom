# -*- coding: utf-8 -*-
#
# Copyright 2022-2023 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Test generation of pleiades cone search used in the paper
"""

import pooch

from rgbloom.__main__ import exec_rgbloom


class Args:
    ra_center = 56.66
    dec_center = 24.10
    search_radius = 1.0
    g_limit = 12.0
    basename = 'rgbloom'
    brightlimit = 8.0
    symbsize = 1.0
    nonumbers = False
    noplot = True            # do not generate PDF file
    nocolor = False
    verbose = False


# compute md5 hash from terminal using:
# linux $ md5sum <filename>
# macOS $ md5 <filename>
auxhash = {
    '200m': "md5:476d8548fbc78f70c7c18806f079cb60",
    'no200m': "md5:b1d6ec03a90824facd1e117a608cd51d"
}

fref = dict()
print('Using reference files:')
for ftype in auxhash:
    fname = f'{Args.basename}_{ftype}.csv'
    ftmp = pooch.retrieve(
        f"http://nartex.fis.ucm.es/~ncl/rgbphot/gaiaDR3/{fname}",
        known_hash=auxhash[ftype]
    )
    fref[ftype] = ftmp
    print(ftmp)


def test_pleiades_simple():
    exec_rgbloom(Args)
    for ftype in auxhash:
        fname = f'{Args.basename}_{ftype}.csv'
        print(f'Checking {fname}')
        with open(fname, 'r') as f:
            file1 = f.read()
        with open(fref[ftype], 'r') as f:
            file2 = f.read()
        assert file1 == file2
