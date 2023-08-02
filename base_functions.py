#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 23:51:30 2023

@author: E0463430
"""

import glob
import gzip
import subprocess
import os
import sys
import shutil
import pandas as pd


########################## TESTING STARTS HERE ############################

def run_fastqc(ARGDICT, dirra_name):
    '''FASTQC on nascnet FASTQ files'''
    if not os.path.isdir('/data/fastqc_output'):
        os.mkdir('/data/fastqc_output')
    if not os.path.isdir(f'/data/fastqc_output/{dirra_name}'):
        os.mkdir(f'/data/fastqc_output/{dirra_name}')
    fastqs_path = glob.glob(f'/data/input_fastq/{dirra_name}/*.fastq.gz')
    fastqs_relevant_path = sorted([i for i in fastqs_path if "_R" in i])
    fastqs_command = " ".join(fastqs_relevant_path)
    command = f'fastqc {fastqs_command} --outdir /data/fastqc_output/{dirra_name} --threads {ARGDICT["nthreads"]}'
    subprocess.call(command, shell=True)


def run_crfastq(ARGDICT):
    '''run cellranger mkfastq'''
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
            run_fastqc(ARGDICT, dirra_name)
    os.chdir('/data')


def run_bclconvert(ARGDICT):
    '''run BCL convert on BCL files'''
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
                    convert_sheet(f'{dirra}/SampleSheet.csv',  ARGDICT["converter_sheet"], nlanes)
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
            run_fastqc(ARGDICT, dirra_name)

                
def organize_crcount(ARGDICT, samples):
    '''organize crcount output'''
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
            in_sample_res = pd.read_csv(
                f'{indira}/outs/metrics_summary.csv', header=0)
            in_sample_res.index = [sample]
            in_sampl = pd.concat([in_sampl, in_sample_res])
    if os.path.isfile('/data/cr_count_organized_output/metrics.csv'):
        outfile = pd.concat([pd.read_csv('/data/cr_count_organized_output/metrics.csv',
                                                   header=0, index_col=0), in_sampl])    
        outfile.to_csv('/data/cr_count_organized_output/metrics.csv')
    else:
        in_sampl.index.name = 'Sample'
        in_sampl.to_csv('/data/cr_count_organized_output/metrics.csv')


def run_crcount(ARGDICT, PACKAGES):
    '''runs cellranger count'''
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
     
    print("ORGANIZING CR COUNT OUTPUT")            
    organize_crcount(ARGDICT, samples)
        

def starsolo_mkref(ARGDICT):
    '''makes a reference genome for starsolo'''
    print('make a reference genome')
    if not os.path.isfile(f'{ARGDICT["star_solo_reference_template"]}/genome.fa'):
        sys.exit('please include genome.fa in the star_solo_genome_generate directory!')
    if not os.path.isfile(f'{ARGDICT["star_solo_reference_template"]}/genes.gtf'):
        sys.exit('please include genes.gtf in the star_solo_reference_template directory!')
    command = (f'STAR --runMode genomeGenerate ' 
               f'--runThreadN {ARGDICT["nthreads"]} '
               f'--genomeDir /data/star_solo_reference '
               f'--genomeFastaFiles {ARGDICT["star_solo_reference_template"]}/genome.fa '
               f'--sjdbGTFfile {ARGDICT["star_solo_reference_template"]}/genes.gtf ')
    subprocess.call(command, shell=True)
    print("REFERENCE GENOME FOR STAR SOLO CREATED SUCCESSFULLY")


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


def organize_starsolo_output(ARGDICT, samples):
    '''organizes STARSolo output into readable format'''
        
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
    

def run_starsolo(ARGDICT, PACKAGES):
    '''this function runs STARSolo'''

    if not os.path.isdir('/data/STAR_output'):
        os.mkdir('/data/STAR_output')
        
    os.chdir('/data/STAR_output')
    
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
            --runThreadN {ARGDICT["nthreads"]} \
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
    organize_starsolo_output(ARGDICT, samples)
