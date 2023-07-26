#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 23:51:30 2023

@author: E0463430
"""

import argparse
import glob
import gzip
import multiprocessing
import subprocess
import os
import sys
import shutil
import pandas as pd


#number of thresds
nthreads = multiprocessing.cpu_count()-1

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
        converter_sheet = get_conversion_tables(ARGDICT["sheet_conversion"])
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

# bcl2fastq: cellranger mkfastq
if "cr_mkfastq" in ARGDICT:
    print("RUNNING CELLRANGER (CR) MKFASTQ")
    if not os.path.isdir('/data/input_bcl'):
        print('please place your bcls in folders in directory input_bcl if you want to convert to fastq')
        sys.exit()
    dirs_to_convert = glob.glob('/data/input_bcl/*/')
    if not os.path.isdir('/data/bcl_conversion'):
        os.mkdir('/data/bcl_conversion')
    if not os.path.isdir('/data/input_fastq'):
        os.mkdir('/data/input_fastq')
    if not dirs_to_convert:
        sys.exit('no directories to convert bcl to fastq')
    for dirra in dirs_to_convert:
        dirra_name = dirra.split('/')[-2]
        print(dirra_name)
        os.chdir('/data/bcl_conversion')
        bcl2fastq_command = ('cellranger mkfastq '
                       f'--run {dirra} '
                       f'--samplesheet {dirra}/SampleSheet.csv '
                       f'--id {dirra_name}')
        print(bcl2fastq_command)
        subprocess.call(f"{bcl2fastq_command}", shell=True)
        if not os.path.isdir(f'/data/input_fastq/{dirra_name}'):
            os.mkdir(f'/data/input_fastq/{dirra_name}')
        outfolders = glob.glob(f'/data/bcl_conversion/{dirra_name}/outs/fastq_path/*/')
        print("iutfolders", outfolders)
        flowcell_folder = [i for i in outfolders if i.split('/')[-2] != "Stats" and i.split('/')[-2] != 'Reports'][0]
        print(flowcell_folder)
        fastq_output = glob.glob(f'{flowcell_folder}/*.fastq.gz')
        print(fastq_output)
        for fastq in fastq_output:
            shutil.move(f'{fastq}', f'/data/input_fastq/{dirra_name}')
        #FASTqc 
        if "fastqc" in ARGDICT:
            print("RUNNING FASTQC")
            if not os.path.isdir('/data/fastqc_output'):
                os.mkdir('/data/fastqc_output')
            if not os.path.isdir(f'/data/fastqc_output/{dirra_name}'):
                os.mkdir(f'/data/fastqc_output/{dirra_name}')
            fastqs_path = glob.glob(f'/data/input_fastq/{dirra_name}/*.fastq.gz')
            fastqs_relevant_path = sorted([i for i in fastqs_path if "_R" in i])
            fastqs_command = " ".join(fastqs_relevant_path)
            command = f'fastqc {fastqs_command} --outdir /data/fastqc_output/{dirra_name} --threads {nthreads}'
            subprocess.call(command, shell=True)
        os.chdir('/data')

if "bcl_convert" in ARGDICT:
    print("RUNNING BCL CONVERT")
    if not os.path.isdir('/data/input_bcl'):
        print('please place your bcls in folders in directory input_bcl if you want to bcl convert')
        sys.exit()
    dirs_to_convert = glob.glob('/data/input_bcl/*/')
    if not os.path.isdir('/data/bcl_conversion'):
        os.mkdir('/data/bcl_conversion')
    if not os.path.isdir('/data/input_fastq'):
        os.mkdir('/data/input_fastq')
    if not dirs_to_convert:
        sys.exit('no directories to convert bcl to fastq')
    for dirra in dirs_to_convert:
        dirra_name = dirra.split('/')[-2]
        # convert sample sheet if needed
        if "bcl_convert_sheet_conv" in ARGDICT:
            from convert_sample_sheet import convert_sheet
            with open(f'{dirra}/SampleSheet.csv') as in_sheet:
                header = in_sheet.readlines()[0]
                nlanes = len(glob.glob(f'{dirra}/Data/Intensities/BaseCalls/*/'))
                if "Lane" in header:
                    convert_sheet(f'{dirra}/SampleSheet.csv', converter_sheet, nlanes)
                    print("SAMPLE SHEET CONVERTED TO BCL CONVERT FORMAT. PLEASE DO NOT CONVERT AGAIN")
                else:
                    print("THE SHEET APPEARS TO NOT BE IN 10X FORMAT. CONVERSION SKIPPED")
        bcl_convert_command = ('bcl-convert '
                       f'--bcl-input-directory {dirra} '
                       f'--output-directory /data/bcl_conversion/{dirra_name} '
                       f'--sample-sheet {dirra}/SampleSheet.csv')
        print(bcl_convert_command)
        subprocess.call(f"{bcl_convert_command}", shell=True)
        if not os.path.isdir(f'/data/input_fastq/{dirra_name}'):
            os.mkdir(f'/data/input_fastq/{dirra_name}')
        else:
            print(f'The fastq directory {dirra_name} already exists. Overriding older files')
        fastq_output = glob.glob(f'/data/bcl_conversion/{dirra_name}/*.fastq.gz')
        for fastq in fastq_output:
            if "Undetermined" not in fastq:
                shutil.move(f'{fastq}', f'/data/input_fastq/{dirra_name}')
        #FASTqc 
        if "fastqc" in ARGDICT:
            print("RUNNING FASTQC")
            if not os.path.isdir('/data/fastqc_output'):
                os.mkdir('/data/fastqc_output')
            if not os.path.isdir(f'/data/fastqc_output/{dirra_name}'):
                os.mkdir(f'/data/fastqc_output/{dirra_name}')
            fastqs_path = glob.glob(f'/data/input_fastq/{dirra_name}/*.fastq.gz')
            fastqs_relevant_path = sorted([i for i in fastqs_path if "_R" in i])
            fastqs_command = " ".join(fastqs_relevant_path)
            command = f'fastqc {fastqs_command} --outdir /data/fastqc_output/{dirra_name} --threads {nthreads}'
            subprocess.call(command, shell=True)
  
# redefine input after potentially converting bcl to fastq
if os.path.isdir('/data/input_fastq/'):
    check_subfolders = glob.glob('/data/input_fastq/*/')
    if len(check_subfolders) > 0:
        ARGDICT["input"] += check_subfolders                

# actually run cellranger
if "cr_count" in ARGDICT or "cr_count_feature" in ARGDICT:
    print("RUNNING CELLRANGER COUNT (CR)")
    # cellranger count
        
    if not os.path.isdir('/data/cr_count_output'):
        os.mkdir('/data/cr_count_output')
        
    os.chdir('/data/cr_count_output')
    
    if "cr_count_feature" not in ARGDICT: # regular count, no feature
        
        if ARGDICT["samp"] == 'All':
            all_fastqs = [glob.glob(f'{locale}/*.fastq.gz')
                          for locale in ARGDICT["input"]]
            all_fastqs = [item for sublist in all_fastqs for item in sublist]
            samples = sorted(set([i.split('/')[-1].split('_')[0] for i in all_fastqs]))
        else:
            samples = ARGDICT["samp"]
        
        all_fastqs_run = ','.join(ARGDICT["input"])
        
        print('samples', samples)
        print('all fastqs run', all_fastqs_run)

        for sample in samples:
            command = (f'{PACKAGES["cellranger"] }/cellranger count '
                       f'--id {sample} '
                       f'--sample {sample} '
                       f'--fastqs {all_fastqs_run} '
                       f'--transcriptome {ARGDICT["cr_count_reference"]} '
                       f'--no-bam ')
            if "cr_count_chemistry" in ARGDICT:
                command += (f'--chemistry {ARGDICT["cr_count_chemistry"]} ')
            if "cr_count_forcecells" in ARGDICT:
                command += (f'--force cells {ARGDICT["cr_count_forcecells"]} ')        
            print(command)
            subprocess.call(command, shell=True)
            
    else:
        if not os.path.isdir('/data/library_files'):
            sys.exit('please put library files names samplecrispr_samplegex.csv inside the library_files directory')
        if not os.path.isfile('/data/feature_ref.csv'):
            sys.exit('please include feature descriptions in feature_ref.csv')
        libraries = sorted(glob.glob('/data/library_files/*.csv'))
        
        samples = []
        for lib_file in libraries:
            sample = lib_file.split('/')[-1].split('.')[0]
            command =  (f'{PACKAGES["cellranger"] }/cellranger count '
                        f'--id {sample} '
                        f'--libraries {lib_file} '
                        f'--feature-ref /data/feature_ref.csv '
                        f'--transcriptome {ARGDICT["cr_count_reference"]} '
                        f'--no-bam')
            if ARGDICT["samp"] == "All" or (ARGDICT["samp"] != "All" and sample in ARGDICT["samp"]):
                samples.append(sample)
                print(command)
                subprocess.call(command, shell=True)
            
    # COMBINE OUTPUT INTO READABLE FORMAT INCL CELLBRIDGE INPUT
    print("ORGANIZING CR OUTPUT")
    if not os.path.isdir('/data/cr_count_organized_output'):
        os.mkdir('/data/cr_count_organized_output')
    if not os.path.isdir('/data/cr_count_organized_output/cellbridge_input'):
        os.mkdir('/data/cr_count_organized_output/cellbridge_input')
    if not os.path.isdir('/data/cr_count_organized_output/loupe_files'):
        os.mkdir('/data/cr_count_organized_output/loupe_files')
    if not os.path.isdir('/data/cr_count_organized_output/web_summaries'):
        os.mkdir('/data/cr_count_organized_output/web_summaries')
    indirs_cellranger = [i for i in glob.glob('/data/cr_count_output/*/') if i.split('/')[-2] in samples]
    in_sampl = pd.DataFrame()
    final_samples = []
    for indira in sorted(indirs_cellranger):
        sample = indira.split('/')[-2]
        if os.path.isfile(f'{indira}/outs/metrics_summary.csv'):
            final_samples.append(sample)
            shutil.copy2(f'{indira}/outs/web_summary.html',
                         f'/data/cr_count_organized_output/web_summaries/{sample}_web_summary.html')
            shutil.copy2(f'{indira}/outs/cloupe.cloupe',
                         f'/data/cr_count_organized_output/loupe_files/{sample}_cloupe.cloupe')
            if not os.path.isdir(f'/data/cr_count_organized_output/cellbridge_input/{sample}'):
                os.mkdir(f'/data/cr_count_organized_output/cellbridge_input/{sample}')
            shutil.copy2(f'{indira}/outs/filtered_feature_bc_matrix.h5',
                         f'/data/cr_count_organized_output/cellbridge_input/{sample}/filtered_feature_bc_matrix.h5')
            in_sampl = pd.concat([in_sampl, pd.read_csv(
                f'{indira}/outs/metrics_summary.csv', header=0)])
    in_sampl.index = final_samples
    in_sampl.index.name = 'Sample'
    in_sampl.to_csv('/data/cr_count_organized_output/metrics.csv')

#STARSolo genome generate
if "star_solo_genome_generate" in ARGDICT:
    print('make a reference genome')
    if not os.path.isfile(f'{ARGDICT["star_solo_reference_template"]}/genome.fa'):
        sys.exit('please include genome.fa in the star_solo_genome_generate directory!')
    if not os.path.isfile(f'{ARGDICT["star_solo_reference_template"]}/genes.gtf'):
        sys.exit('please include genes.gtf in the star_solo_reference_template directory!')
    command = (f'STAR --runMode genomeGenerate ' 
               f'--runThreadN {nthreads} '
               f'--genomeDir /data/star_solo_reference '
               f'--genomeFastaFiles {ARGDICT["star_solo_reference_template"]}/genome.fa '
               f'--sjdbGTFfile {ARGDICT["star_solo_reference_template"]}/genes.gtf ')
    subprocess.call(command, shell=True)
    print("REFERENCE GENOME FOR STAR SOLO CREATED SUCCESSFULLY")

# define the reference dirs
if os.path.isdir("/data/star_solo_reference"):    
    ARGDICT["star_solo_reference"] = "/data/star_solo_reference"

# quick check
if "star_solo" in ARGDICT and "star_solo_reference" not in ARGDICT:
    sys.exit("please provide a reference genome \
             for STAR Solo in the star_solo_reference directory")

#STARSolo
if "star_solo" in ARGDICT:
    print("RUNNING STAR SOLO")
                
    if not os.path.isdir('/data/STAR_output'):
        os.mkdir('/data/STAR_output')
        
    os.chdir('/data/STAR_output')
        
    def get_fastq_info(infiles):
        '''goes over fastq to prepare star solo commands'''
        infiles_sorted = {}
        for filey in infiles:
            filename = filey.split('/')[-1]
            name = filename.split('_')[0]
            if name not in infiles_sorted:
                infiles_sorted[name] = {}
                infiles_sorted[name]['R2'] = []
                infiles_sorted[name]['R1'] = []
            if 'R1' in filey:
                infiles_sorted[name]['R1'].append(filey)
            if 'R2' in filey:
                infiles_sorted[name]['R2'].append(filey)
        for namey in infiles_sorted:
            infiles_sorted[namey]['R1'] = ",".join(sorted(infiles_sorted[namey]['R1']))
            infiles_sorted[namey]['R2'] = ",".join(sorted(infiles_sorted[namey]['R2']))
        return infiles_sorted
    
    infiles = []
    for indir in ARGDICT["input"]:
        infiles += glob.glob(f'{indir}/*.fastq.gz')
    
    # get info for all the samples
    # if samples wre to be excluded, will not iterate over them (below)
    starsolo_fastq_info = get_fastq_info(infiles)
    
    if ARGDICT["samp"] == 'All':
         samples = sorted(starsolo_fastq_info.keys())
    else:
         samples = ARGDICT["samp"]
    
    for namey in samples:
        staroutputdir = f'/data/STAR_output/{namey}'        

        if ARGDICT["star_solo_chem"] not in ARGDICT["chemistry_file"].index:
            sys.exit(f'The specified chemistry {ARGDICT["star_solo_chem"]} is not available \
                     Please specify one of {",".join(ARGDICT["chemistry_file"].index.to_list())}')
        
        barcode_file = ARGDICT["chemistry_file"].loc[ARGDICT["star_solo_chem"], "List"]
        barcode_len = ARGDICT["chemistry_file"].loc[ARGDICT["star_solo_chem"], "barcode"]
        umi_len = ARGDICT["chemistry_file"].loc[ARGDICT["star_solo_chem"],"UMI"]
        
        if namey in starsolo_fastq_info:
        
            command = (f'{PACKAGES["star_path"]}/STAR \
            --runThreadN {nthreads} \
            --genomeDir {ARGDICT["star_solo_reference"]} \
            --readFilesIn {starsolo_fastq_info[namey]["R2"]} {starsolo_fastq_info[namey]["R1"]} \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist {ARGDICT["barcodes"]}/{barcode_file} \
            --readFilesCommand zcat \
            --soloCBlen {barcode_len} \
            --soloUMIlen {umi_len} \
            --outFileNamePrefix {staroutputdir} \
            --soloFeatures {ARGDICT["star_solo_features"]}')
            
            print(command)
            subprocess.call(command, shell=True)
        
        else:
            print(f'Sample {namey} is not available! Please check the name and retry')      
        
    print("ORGANIZING SOLO OUTPUT")
    if not os.path.isdir('/data/STAR_organized_output'):
        os.mkdir('/data/STAR_organized_output')
        
    if not os.path.isdir('/data/STAR_organized_output/cellbridge_input'):
        os.mkdir('/data/STAR_organized_output/cellbridge_input')
        
    indirs_star = [i for i in glob.glob('/data/STAR_output/*/') if i.split('/')[-2].split('Solo.out')[0] in samples]
    
    in_sampl = pd.DataFrame()
    
    for indira in sorted(indirs_star):
        sample_old = indira.split('/')[-2]
        sample_new = sample_old.split('Solo.out')[0]
        print(f'processing {sample_new}')
        if os.path.isfile(f'{indira}/{ARGDICT["star_solo_features"]}/Summary.csv'):
            
            if not os.path.isdir(f'/data/STAR_organized_output/cellbridge_input/{sample_new}'):
                os.mkdir(f'/data/STAR_organized_output/cellbridge_input/{sample_new}')
                
            in_mtx = f'{indira}/{ARGDICT["star_solo_features"]}/filtered/matrix.mtx'
            out_mtx = f'/data/STAR_organized_output/cellbridge_input/{sample_new}/matrix.mtx.gz'
            with open(in_mtx, 'rb') as f_in, gzip.open(out_mtx, 'wb', compresslevel=6) as f_out:
                  f_out.write(f_in.read())
            in_feat = f'{indira}/{ARGDICT["star_solo_features"]}/filtered/features.tsv'
            out_feat = f'/data/STAR_organized_output/cellbridge_input/{sample_new}/features.tsv.gz'
            with open(in_feat, 'rb') as f_in, gzip.open(out_feat, 'wb', compresslevel=6) as f_out:
                  f_out.write(f_in.read())
            in_bar = f'{indira}/{ARGDICT["star_solo_features"]}/filtered/barcodes.tsv'
            out_bar = f'/data/STAR_organized_output/cellbridge_input/{sample_new}/barcodes.tsv.gz'
            with open(in_bar, 'rb') as f_in, gzip.open(out_bar, 'wb', compresslevel=6) as f_out:
                  f_out.write(f_in.read())
                  
            overview = pd.read_csv(f'{indira}/{ARGDICT["star_solo_features"]}/Summary.csv',
                               index_col=0, header=None)
            overview.columns = [sample_new]
            overview = overview.T
            
            in_sampl = pd.concat([in_sampl, overview])
      
    in_sampl.index.name = 'Sample'
    in_sampl.to_csv('/data/STAR_organized_output/metrics.csv')
  
 #########################################################################