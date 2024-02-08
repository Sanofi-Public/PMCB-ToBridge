#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 23:17:45 2023

@author: E0463430
"""

import pandas as pd
import glob
import shutil

def get_conversion_tables(path_to_dir):
    '''gets all the conversion tables for the sample sheet conversion'''
    all_conv = []
    conversion_csvs = glob.glob(f'{path_to_dir}/*.csv')
    for conv in conversion_csvs:
        convert = pd.read_csv(conv, index_col=0, header=3)
        all_conv.append(convert)
    combino = pd.concat(all_conv)
    return combino

# this needs to account for a different number of lanes
# will deal with this
def convert_sheet(infile, combino, nlanes):
    '''converts sample sheet from 10x format to bcl convert friendly format'''
    
    outfile = infile # override the original file
    shutil.copy2(infile, "/".join(infile.split("/")[:-1]) + "/SampleSheet10x.csv")
    infile = "/".join(infile.split("/")[:-1]) + "/SampleSheet10x.csv"

    header = ('[Header]\n'
              'FileFormatVersion,2\n\n'
              '[BCLConvert_Settings]\n'
              'CreateFastqForIndexReads,0\n\n'
              '[BCLConvert_Data]\n'
              'Lane,Sample_ID,index,index2\n')

    with open(outfile, "w") as out_handle:
        out_handle.write(header)
        with open(infile) as in_handle:
            for liney in in_handle:
                parsed_liney = liney.rstrip().split(',')
                if 'Lane' not in parsed_liney:
                    if parsed_liney[2] in combino.index:
                        for lane in range(1, nlanes+1):
                            out_handle.write(f'{lane},{parsed_liney[1]},{combino.loc[parsed_liney[2], "index(i7)"]},{combino.loc[parsed_liney[2], "index2_workflow_b(i5)"]}'+'\n')
                    else:
                        print(f'warning! index {parsed_liney[2]} cannot be converted. please make sure you use the correct and most up-to-date csv files')
            