#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 23:51:30 2023

@author: E0463430
"""

import argparse
import glob
import multiprocessing
import os
import sys
import pandas as pd
import base_functions as bf


# ARGUMENTS
PARSER = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter)

PARSER.add_argument("--cr_mkfastq", required=False, action='store_true', default=False,
                    help="Run cellranger make fastq")

PARSER.add_argument("--bcl_convert", required=False, action='store_true', default=False,
                    help="Run BCL convert")
PARSER.add_argument("--bcl_convert_sheet_conv", required=False, action='store_true', default=False,
                    help="Convert sample sheet to use BCL convert")

PARSER.add_argument("--fastqc", required=False, action='store_true', default=True,
                    help="Run FASTQC")

PARSER.add_argument("--star_solo", required=False, action='store_true', default=False,
                    help="Path to STAR Solo input")
PARSER.add_argument("--star_solo_chem", required=False, 
                    help="What is the chemistry?")
PARSER.add_argument("--star_solo_features", required=False, default="GeneFull_Ex50pAS",
                    help="Quantification of different features")
PARSER.add_argument("--star_solo_genome_generate", action='store_true', required=False, 
                    help="Generate a ref genome")

PARSER.add_argument("--cr_count", required=False, action='store_true', default=False,
                    help="Path to cellranger input")
PARSER.add_argument("--cr_count_chemistry", required=False,
                    help="Specify chemistry if cannot be detected automatically")
PARSER.add_argument("--cr_count_forcecells", required=False,
                    help="Force a certain number of cells")

PARSER.add_argument("--cr_count_feature", required=False, action='store_true', default=False,
                    help="Feature bacrode assay")

PARSER.add_argument("--samp", required=False, default="All",
                    help="If you only want cellranger run on specific samples, use this flag")

PARSER.add_argument("--nthreads", required=False, default=multiprocessing.cpu_count()-1,
                    help="Number of threads")

# PACKAGES
PACKAGES = {}

PACKAGES['bcl2fastq'] = '/usr/local/bin/bcl2fastq'
os.environ['PATH'] += os.pathsep + os.pathsep.join([PACKAGES['bcl2fastq']])

PACKAGES['bcl_convert'] = '/usr/local/bin/bcl-convert'
os.environ['PATH'] += os.pathsep + os.pathsep.join([PACKAGES['bcl_convert']])

PACKAGES['cellranger'] = '/opt/tobridge/cellranger-7.1.0/bin'
os.environ['PATH'] += os.pathsep + os.pathsep.join([PACKAGES['cellranger']])

PACKAGES['fastq_path'] = '/opt/tobridge/FastQC'
os.environ['PATH'] += os.pathsep + os.pathsep.join([PACKAGES['fastq_path']])

PACKAGES['star_path'] = '/opt/tobridge/STAR-2.7.10b/source'
os.environ['PATH'] += os.pathsep + os.pathsep.join([PACKAGES['star_path']])

# HARD-CODED PARAMETERS FIRST
# nota that STAR Solo reference is coded in below

ARGDICT = {}
ARGDICT["input"] = []

if os.path.isdir("/data/cr_count_reference"):
    ARGDICT["cr_count_reference"] = "/data/cr_count_reference"
    
if os.path.isdir("/data/star_solo_reference_template"):    
    ARGDICT["star_solo_reference_template"] = "/data/star_solo_reference_template"

ARGDICT["barcodes"] = "/opt/tobridge/barcodes"
ARGDICT["sheet_conversion"] = "/opt/tobridge/index_csvs"
    
ARGDICT["chemistry_file"] = pd.read_csv("/opt/tobridge/chemistry.csv",
                                        header=0,
                                        index_col=0)

# PASSED-IN ARGUMENTS
ARGIES = PARSER.parse_args()

if ARGIES.cr_mkfastq:
    ARGDICT["cr_mkfastq"] = ARGIES.cr_mkfastq

if ARGIES.bcl_convert:
    ARGDICT["bcl_convert"] = ARGIES.bcl_convert
    
if "bcl_convert" in ARGDICT and "cr_mkfastq" in ARGDICT:
    sys.exit("please choose either bcl convert or cellranger mkfastq, but not both")
    
if ARGIES.bcl_convert_sheet_conv:
    ARGDICT["bcl_convert_sheet_conv"] = ARGIES.bcl_convert_sheet_conv
    from convert_sample_sheet import get_conversion_tables
    if os.path.isdir(ARGDICT["sheet_conversion"]):
        ARGDICT["converter_sheet"] = get_conversion_tables(ARGDICT["sheet_conversion"])
    else:
        sys.exit("nothing to convert ")
    
if ARGIES.fastqc:
    ARGDICT["fastqc"] = ARGIES.fastqc

if ARGIES.cr_count:
    ARGDICT["cr_count"] = ARGIES.cr_count
    
if ARGIES.cr_count_chemistry:
    ARGDICT["cr_count_chemistry"] = ARGIES.cr_count_chemistry

if ARGIES.cr_count_forcecells:
    ARGDICT["cr_count_forcecells"] = ARGIES.cr_count_forcecells
    
if ARGIES.cr_count_feature:
    ARGDICT["cr_count_feature"] = ARGIES.cr_count_feature

if ("cr_count" in ARGDICT or "cr_count_feature" in ARGDICT) and "cr_count_reference" not in ARGDICT:
    sys.exit("please provide a reference genome \
             for cellranger count in the cr_count_reference directory")
    
if ARGIES.star_solo:
    ARGDICT["star_solo"] = ARGIES.star_solo
    
if ARGIES.star_solo_genome_generate:
    ARGDICT["star_solo_genome_generate"] = ARGIES.star_solo_genome_generate
    if "star_solo_reference_template" not in ARGDICT:
        sys.exit("Directory star_solo_reference_template "
                 "containing genome.fa and genes.gtf "
                 "is NOT present. Cannot make reference for STARSolo.")
    
if "star_solo" in ARGDICT and not ARGIES.star_solo_chem:
    sys.exit('Please provide chemistry for your STAR Solo run')
    
if ARGIES.star_solo_chem:
    ARGDICT["star_solo_chem"] = ARGIES.star_solo_chem

ARGDICT["star_solo_features"] = ARGIES.star_solo_features

if ARGIES.samp == 'All':
    ARGDICT["samp"] = 'All'
else:
    ARGDICT["samp"] = ARGIES.samp.split(',')
    
ARGDICT["nthreads"] = ARGIES.nthreads

######################### RUN SCRIPT ##########################

if __name__ == "__main__":
    
    if "cr_mkfastq" in ARGDICT:
        print("RUNNING CELLRANGER (CR) MKFASTQ")
        bf.run_crfastq(ARGDICT)
    
    if "bcl_convert" in ARGDICT:
        print("RUNNING BCL CONVERT")
        bf.run_bclconvert(ARGDICT)
      
    # redefine input after potentially converting bcl to fastq
    if os.path.isdir('/data/input_fastq/'):
        check_subfolders = glob.glob('/data/input_fastq/*/')
        if len(check_subfolders) > 0:
            ARGDICT["input"] += check_subfolders     
    
    # actually run cellranger
    if "cr_count" in ARGDICT or "cr_count_feature" in ARGDICT:
        print("RUNNING CELLRANGER COUNT (CR) AND ORGANIZING OUTPUT")
        bf.run_crcount(ARGDICT, PACKAGES)
     
    #STARSolo genome generate
    if "star_solo_genome_generate" in ARGDICT:
        bf.starsolo_mkref(ARGDICT)
    
    # define the reference dirs
    if os.path.isdir("/data/star_solo_reference"):    
        ARGDICT["star_solo_reference"] = "/data/star_solo_reference"
    
    # quick check
    if "star_solo" in ARGDICT and "star_solo_reference" not in ARGDICT:
        sys.exit("please provide a reference genome \
                 for STAR Solo in the star_solo_reference directory")
    
    #STARSolo
    if "star_solo" in ARGDICT:
        print("RUNNING STAR SOLO AND ORGANIZING OUTPUT")
        bf.run_starsolo(ARGDICT, PACKAGES)
      
 #############################################################################