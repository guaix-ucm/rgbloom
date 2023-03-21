# -*- coding: utf-8 -*-
#
# Copyright 2022-2023 Universidad Complutense de Madrid
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Generate output CSV files
"""


def step5(r_dr3_200m, r_dr3_no200m, basename, verbose):
    """Generate output CSV files

    Parameters
    ----------
    r_dr3_200m : astropy Table
        Objects in both the DR3 query and the 200M sample.
    r_dr3_no200m : astropy Table
        Objects in the DR3 query no present in the 200M sample.
    basename : str
        Base name for output files.
    verbose : bool
        If True, display additional information.

    """

    print('<STEP5> Saving output CSV files')

    # ---
    # RGB values from r_dr3_200m

    # columns to be saved (use a list to guarantee the same order)
    outcolumns_list = [
        'source_id', 'ra', 'dec',
        'RGB_B', 'RGB_G', 'RGB_R',
        'errRGB_B', 'errRGB_G', 'errRGB_R',
        'objtype', 'qlflag'
    ]
    # define column format with a dictionary
    outcolumns = {
        'source_id': '19d',
        'ra': '14.9f',
        'dec': '14.9f',
        'RGB_B': '8.4f',
        'RGB_G': '8.4f',
        'RGB_R': '8.4f',
        'errRGB_B': '8.4f',
        'errRGB_G': '8.4f',
        'errRGB_R': '8.4f',
        'objtype': '1d',
        'qlflag': '1d'
    }
    if set(outcolumns_list) != set(outcolumns.keys()):
        raise SystemExit('ERROR: check outcolumns_list and outcolumns')
    csv_header = 'number,' + ','.join(outcolumns_list)
    outfile = f'{basename}_200m.csv'
    f = open(outfile, 'wt')
    f.write(csv_header + '\n')
    for row in r_dr3_200m:
        cout = []
        for item in outcolumns_list:
            cout.append(eval("f'{row[item]:" + f'{outcolumns[item]}' + "}'"))
        f.write(f"{row['number']:8d}, " + ','.join(cout) + '\n')
    f.close()
    print(f'        --> file {outfile} saved')

    # ---
    # RGB values from r_dr3_no200m

    # columns to be saved (use a list to guarantee the same order)
    outcolumns_list = [
        'source_id', 'ra', 'dec',
        'phot_variable_flag',
        'bp_rp',
        'RGB_B_C21', 'RGB_G_C21', 'RGB_R_C21'
    ]
    # define column format with a dictionary
    outcolumns = {
        'source_id': '19d',
        'ra': '14.9f',
        'dec': '14.9f',
        'phot_variable_flag': '>15',
        'bp_rp': '7.3f',
        'RGB_B_C21': '7.3f',
        'RGB_G_C21': '7.3f',
        'RGB_R_C21': '7.3f'
    }
    if set(outcolumns_list) != set(outcolumns.keys()):
        raise SystemExit('ERROR: check outcolumns_list and outcolumns')
    csv_header = 'number,' + ','.join(outcolumns_list)
    outfile = f'{basename}_no200m.csv'
    f = open(outfile, 'wt')
    f.write(csv_header + '\n')
    for row in r_dr3_no200m:
        cout = []
        for item in outcolumns_list:
            cout.append(eval("f'{row[item]:" + f'{outcolumns[item]}' + "}'"))
        f.write(f"{row['number']:8d}, " + ','.join(cout) + '\n')
    f.close()
    print(f'        --> file {outfile} saved')
