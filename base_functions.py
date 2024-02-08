#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 23:51:30 2023

@author: Andre Kurlovs and Nima Nouri
"""

import csv
import glob
import gzip
import subprocess
import os
import sys
import shutil
import pandas as pd
from multiprocessing import Pool
import pysam
from datetime import datetime
import psutil

########################## TESTING STARTS HERE ############################

def run_fastqc(ARGDICT):
    '''FASTQC on nascnet FASTQ files'''
    if not os.path.isdir('/data/fastqc_output'):
        os.mkdir('/data/fastqc_output')
    for dirra in glob.glob('/data/input_fastq/*/'):
        dirra_name = dirra.split('/')[-2]
        if not os.path.isdir(f'/data/fastqc_output/{dirra_name}'):
            os.mkdir(f'/data/fastqc_output/{dirra_name}')
        fastqs_path = glob.glob(f'/data/input_fastq/{dirra_name}/*.fastq.gz')
        fastqs_relevant_path = sorted([i for i in fastqs_path if "_R" in i])
        fastqs_command = " ".join(fastqs_relevant_path)
        command = f'fastqc {fastqs_command} --outdir /data/fastqc_output/{dirra_name} --threads {ARGDICT["nthreads"]}'
        subprocess.run(command, shell=True)


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
        subprocess.run(f"{bcl2fastq_command}", shell=True)
        if not os.path.isdir(f'/data/input_fastq/{dirra_name}'):
            os.mkdir(f'/data/input_fastq/{dirra_name}')
        outfolders = glob.glob(f'/data/bcl_conversion/{dirra_name}/outs/fastq_path/*/')
        print("outfolders", outfolders)
        flowcell_folder = [i for i in outfolders if i.split('/')[-2] != "Stats" and i.split('/')[-2] != 'Reports'][0]
        if len(flowcell_folder) == 0:
            sys.exit('conversion did not take place successfully')
        fastq_output = glob.glob(f'{flowcell_folder}/*.fastq.gz')
        print(fastq_output)
        for fastq in fastq_output:
            shutil.move(f'{fastq}', f'/data/input_fastq/{dirra_name}')
    os.chdir('/data')


def run_bclconvert(ARGDICT):
    '''run BCL convert on BCL files'''
    if not os.path.isdir('/data/input_bcl'):
        print('please place your bcls in folders in directory input_bcl if you want to bcl convert')
        sys.exit()
    dirs_to_convert = glob.glob('/data/input_bcl/*/')
    if not os.path.isdir('/data/bcl_conversion'):
        os.mkdir('/data/bcl_conversion')
    else:
        shutil.rmtree('/data/bcl_conversion')
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
                    convert_sheet(f'{dirra}/SampleSheet.csv', ARGDICT["converter_sheet"], nlanes)
                    print("SAMPLE SHEET CONVERTED TO BCL CONVERT FORMAT. PLEASE DO NOT CONVERT AGAIN")
                else:
                    print("THE SHEET APPEARS TO NOT BE IN 10X FORMAT. CONVERSION SKIPPED")
        bcl_convert_command = ('bcl-convert '
                       f'--bcl-input-directory {dirra} '
                       f'--output-directory /data/bcl_conversion/{dirra_name} '
                       f'--sample-sheet {dirra}/SampleSheet.csv')
        print(bcl_convert_command)
        subprocess.run(f"{bcl_convert_command}", shell=True)
        if not os.path.isdir(f'/data/input_fastq/{dirra_name}'):
            os.mkdir(f'/data/input_fastq/{dirra_name}')
        else:
            print(f'The fastq directory {dirra_name} already exists. Overriding older files')
        fastq_output = glob.glob(f'/data/bcl_conversion/{dirra_name}/*.fastq.gz')
        for fastq in fastq_output:
            if "Undetermined" not in fastq:
                shutil.move(f'{fastq}', f'/data/input_fastq/{dirra_name}')


def parse_fasta(file_path):
    '''parse a FASTA file'''
    sequences = {}
    current_id = None
    current_seq = ""

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = current_seq
                current_id = line[1:]
                current_seq = ""
            else:
                current_seq += line

        if current_id is not None and current_seq:
            sequences[current_id] = current_seq

    return sequences


def run_crmkref(ARGDICT):
    '''make a reference genome for cell ranger'''

    if not os.path.isfile(f'{ARGDICT["cr_count_reference_template"]}/genome.fa'):
        sys.exit('please include genome.fa in the cr_count_reference_template directory!')
    if not os.path.isfile(f'{ARGDICT["cr_count_reference_template"]}/genes.gtf'):
        sys.exit('please include genes.gtf in the cr_count_reference_template directory!')
        
    # Example usage:
    original_gtf = f'{ARGDICT["cr_count_reference_template"]}/genes.gtf'
    gtf_seq = 'GENE\tunknown\texon\t1\tTOTALLEN\t.\t+\t.\tgene_id "GENE"; transcript_id "GENE"; gene_name "GENE"; gene_biotype "protein_coding";\n'
    
    original_fasta = f'{ARGDICT["cr_count_reference_template"]}/genome.fa'  # Replace this with the path to your FASTA file
    additional_fasta = f'{ARGDICT["cr_count_reference_template"]}/add.fa'  # Replace this with the path to your FASTA file

    #unless genome is getting modified, the originals are the finals
    final_fasta = original_fasta
    final_gtf = original_gtf
    
    if os.path.isfile(additional_fasta):
        
        print('incorporating additional sequences')
    
        final_gtf = os.path.join(os.path.dirname(original_gtf), 'genes.gtf.modified')
        final_fasta = os.path.join(os.path.dirname(original_fasta), 'genome.fa.modified')
        
        shutil.copy2(f'{original_gtf}', f'{final_gtf}')
        shutil.copy2(f'{original_fasta}', f'{final_fasta}')
        
        sequences_dict = parse_fasta(additional_fasta)
        
        for gene in sequences_dict:
            gtf_out = gtf_seq.replace('GENE', gene)
            gtf_out = gtf_out.replace('TOTALLEN', f'{len(sequences_dict[gene])}')
            fasta_out = f'>{gene}\n{sequences_dict[gene]}\n'
            with open(final_fasta, "a") as fasta_append:
                fasta_append.write(fasta_out)
            with open(final_gtf, "a") as gtf_append:
                gtf_append.write(gtf_out)
                
    os.chdir('/data/cr_count_reference_template/')
            
    command = (f'cellranger mkref  '
               f'--genome {ARGDICT["cr_genome_generate"]} '
               f'--fasta {final_fasta} '
               f'--genes {final_gtf}')
    
    subprocess.run(command, shell=True)
    shutil.copytree(f'/data/cr_count_reference_template/{ARGDICT["cr_genome_generate"]}',
                     '/data/cr_count_reference')
    os.chdir('/data')
    

def organize_crcount(ARGDICT, samples):
    '''organize crcount output'''
    # COMBINE OUTPUT INTO READABLE FORMAT INCL CELLBRIDGE INPUT
    if not os.path.isdir('/data/cr_count_organized_output'):
        os.mkdir('/data/cr_count_organized_output')
    if not os.path.isdir('/data/cr_count_organized_output/cellbridge'):
        os.mkdir('/data/cr_count_organized_output/cellbridge')
    if not os.path.isdir('/data/cr_count_organized_output/loupe_files'):
        os.mkdir('/data/cr_count_organized_output/loupe_files')
    if not os.path.isdir('/data/cr_count_organized_output/web_summaries'):
        os.mkdir('/data/cr_count_organized_output/web_summaries')
    if 'cr_count_bam' in ARGDICT and not os.path.isdir('/data/cr_count_organized_output/bam_files'):
        os.mkdir('/data/cr_count_organized_output/bam_files')
    indirs_cellranger = [i for i in glob.glob('/data/cr_count_output/*/') if i.split('/')[-2] in samples]
    in_sampl = pd.DataFrame()
    final_samples = []
    for indira in sorted(indirs_cellranger):
        sample = indira.split('/')[-2]
        if os.path.isfile(f'{indira}/outs/metrics_summary.csv'):
            final_samples.append(sample)
            if 'cr_count_bam' in ARGDICT:
                shutil.copy2(f'{indira}/outs/possorted_genome_bam.bam',
                     f'/data/cr_count_organized_output/bam_files/{sample}_possorted_genome_bam.bam')
                shutil.copy2(f'{indira}/outs/possorted_genome_bam.bam.bai',
                     f'/data/cr_count_organized_output/bam_files/{sample}_possorted_genome_bam.bam.bai')
            shutil.copy2(f'{indira}/outs/web_summary.html',
                         f'/data/cr_count_organized_output/web_summaries/{sample}_web_summary.html')
            shutil.copy2(f'{indira}/outs/cloupe.cloupe',
                         f'/data/cr_count_organized_output/loupe_files/{sample}_cloupe.cloupe')
            if not os.path.isdir(f'/data/cr_count_organized_output/cellbridge/{sample}'):
                os.mkdir(f'/data/cr_count_organized_output/cellbridge/{sample}')
            shutil.copy2(f'{indira}/outs/filtered_feature_bc_matrix.h5',
                         f'/data/cr_count_organized_output/cellbridge/{sample}/filtered_feature_bc_matrix.h5')
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


def run_crcount(ARGDICT):
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
            samples = sorted(set([i.split('/')[-1].split('_S')[0] for i in all_fastqs]))
        else:
            samples = ARGDICT["samp"]
        
        all_fastqs_run = ','.join(ARGDICT["input"])
        
        print('samples', samples)
        print('all fastqs run', all_fastqs_run)

        for sample in samples:
            command = (f'cellranger count '
                       f'--id {sample} '
                       f'--sample {sample} '
                       f'--fastqs {all_fastqs_run} '
                       f'--transcriptome {ARGDICT["cr_count_reference"]} '
                       f'--no-bam ')
            if "cr_count_chemistry" in ARGDICT:
                command += (f'--chemistry {ARGDICT["cr_count_chemistry"]} ')
            if "cr_count_forcecells" in ARGDICT:
                command += (f'--force-cells {ARGDICT["cr_count_forcecells"]} ')
            if "cr_count_bam" in ARGDICT:
                command = command.replace(' --no-bam', '')
            print(command)
            subprocess.run(command, shell=True)
            
    else:
        if not os.path.isdir('/data/library_files'):
            sys.exit('please put library files names samplecrispr_samplegex.csv inside the library_files directory')
        if not os.path.isfile('/data/feature_ref.csv'):
            sys.exit('please include feature descriptions in feature_ref.csv')
        libraries = sorted(glob.glob('/data/library_files/*.csv'))
        
        samples = []
        for lib_file in libraries:
            sample = os.path.dirname(lib_file).split('.')[0]
            command =  (f'cellranger count '
                        f'--id {sample} '
                        f'--libraries {lib_file} '
                        f'--feature-ref /data/feature_ref.csv '
                        f'--transcriptome {ARGDICT["cr_count_reference"]} '
                        f'--no-bam')
            if "cr_count_chemistry" in ARGDICT:
                command += (f'--chemistry {ARGDICT["cr_count_chemistry"]} ')
            if "cr_count_forcecells" in ARGDICT:
                command += (f'--force-cells {ARGDICT["cr_count_forcecells"]} ')
            if "cr_count_bam" in ARGDICT:
                command = command.replace(' --no-bam', '')
            if ARGDICT["samp"] == "All" or (ARGDICT["samp"] != "All" and sample in ARGDICT["samp"]):
                samples.append(sample)
                print(command)
                subprocess.run(command, shell=True)
     
    print("ORGANIZING CR COUNT OUTPUT")       
    organize_crcount(ARGDICT, samples)
    

def starsolo_mkref(ARGDICT):
    '''makes a reference genome for starsolo'''
    print('make a reference genome')
    if not os.path.isfile(f'{ARGDICT["star_solo_reference_template"]}/genome.fa'):
        sys.exit('please include genome.fa in the star_solo_genome_generate directory!')
    if not os.path.isfile(f'{ARGDICT["star_solo_reference_template"]}/genes.gtf'):
        sys.exit('please include genes.gtf in the star_solo_reference_template directory!')
        
    # Example usage:
    original_gtf = f'{ARGDICT["star_solo_reference_template"]}/genes.gtf'
    gtf_seq = 'GENE\tunknown\texon\t1\tTOTALLEN\t.\t+\t.\tgene_id "GENE"; transcript_id "GENE"; gene_name "GENE"; gene_biotype "protein_coding";\n'
    
    original_fasta = f'{ARGDICT["star_solo_reference_template"]}/genome.fa'  # Replace this with the path to your FASTA file
    additional_fasta = f'{ARGDICT["star_solo_reference_template"]}/add.fa'  # Replace this with the path to your FASTA file
    
    # unless modification takes place, the finals are the originals
    final_fasta = original_fasta
    final_gtf = original_gtf
    
    if os.path.isfile(additional_fasta):
        
        print('incorporating additional sequences')
    
        final_gtf = os.path.join(os.path.dirname(original_gtf), 'genes.gtf.modified')
        final_fasta = os.path.join(os.path.dirname(original_fasta), 'genome.fa.modified')
        
        shutil.copy2(f'{original_gtf}', f'{final_gtf}')
        shutil.copy2(f'{original_fasta}', f'{final_fasta}')
        
        sequences_dict = parse_fasta(additional_fasta)
        
        for gene in sequences_dict:
            gtf_out = gtf_seq.replace('GENE', gene)
            gtf_out = gtf_out.replace('TOTALLEN', f'{len(sequences_dict[gene])}')
            fasta_out = f'>{gene}\n{sequences_dict[gene]}\n'
            with open(final_fasta, "a") as fasta_append:
                fasta_append.write(fasta_out)
            with open(final_gtf, "a") as gtf_append:
                gtf_append.write(gtf_out)
               
    command = (f'STAR --runMode genomeGenerate ' 
               f'--runThreadN {ARGDICT["nthreads"]} '
               f'--genomeDir /data/star_solo_reference '
               f'--genomeFastaFiles {final_fasta} '
               f'--sjdbGTFfile {final_gtf} ')
    subprocess.run(command, shell=True)
    print("REFERENCE GENOME FOR STAR SOLO CREATED SUCCESSFULLY")


def get_fastq_info(infiles):
    '''goes over fastq to prepare star solo commands'''
    infiles_sorted = {}
    for filey in infiles:
        filename = os.path.basename(filey)
        name = filename.split('_S')[0]
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
    
    if 'SJ' in ARGDICT['star_solo_features']:
        # fix broken symlinks
        for tabfile in glob.glob('/data/STAR_output/*SJ.out.tab'):
            outdir=tabfile.split('SJ.out.tab')[0]+'Solo.out'
            os.remove(f'{outdir}/SJ/raw/features.tsv')
            shutil.copy2(tabfile, f'{outdir}/SJ/raw/features.tsv')
    
    if 'GeneFull_Ex50pAS' not in ARGDICT['star_solo_features']:
        sys.exit('only GeneFull_Ex50pAS gets organization')
        
    if not os.path.isdir('/data/STAR_organized_output'):
        os.mkdir('/data/STAR_organized_output')
        
    if not os.path.isdir('/data/STAR_organized_output/cellbridge'):
        os.mkdir('/data/STAR_organized_output/cellbridge')
        
    indirs_star = [i for i in glob.glob('/data/STAR_output/*/') if i.split('/')[-2].split('Solo.out')[0] in samples]
    
    in_sampl = pd.DataFrame()
    
    for indira in sorted(indirs_star):
        sample_old = indira.split('/')[-2]
        sample_new = sample_old.split('Solo.out')[0]
        print(f'processing {sample_new}')
        if os.path.isfile(f'{indira}/GeneFull_Ex50pAS/Summary.csv'):
            
            if not os.path.isdir(f'/data/STAR_organized_output/cellbridge/{sample_new}'):
                os.mkdir(f'/data/STAR_organized_output/cellbridge/{sample_new}')
                
            in_mtx = f'{indira}/GeneFull_Ex50pAS/filtered/matrix.mtx'
            out_mtx = f'/data/STAR_organized_output/cellbridge/{sample_new}/matrix.mtx.gz'
            with open(in_mtx, 'rb') as f_in, gzip.open(out_mtx, 'wb', compresslevel=6) as f_out:
                  f_out.write(f_in.read())
            in_feat = f'{indira}/GeneFull_Ex50pAS/filtered/features.tsv'
            out_feat = f'/data/STAR_organized_output/cellbridge/{sample_new}/features.tsv.gz'
            with open(in_feat, 'rb') as f_in, gzip.open(out_feat, 'wb', compresslevel=6) as f_out:
                  f_out.write(f_in.read())
            in_bar = f'{indira}/GeneFull_Ex50pAS/filtered/barcodes.tsv'
            out_bar = f'/data/STAR_organized_output/cellbridge/{sample_new}/barcodes.tsv.gz'
            with open(in_bar, 'rb') as f_in, gzip.open(out_bar, 'wb', compresslevel=6) as f_out:
                  f_out.write(f_in.read())
                  
            overview = pd.read_csv(f'{indira}/GeneFull_Ex50pAS/Summary.csv',
                               index_col=0, header=None)
            overview.columns = [sample_new]
            overview = overview.T
            
            in_sampl = pd.concat([in_sampl, overview])

    in_sampl.index.name = 'Sample'
    in_sampl.to_csv('/data/STAR_organized_output/metrics.csv')
    

def run_starsolo(ARGDICT):
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
        
        barcode_file = ARGDICT["chemistry_file"].loc[ARGDICT["star_solo_chem"], "List"]
        barcode_len = ARGDICT["chemistry_file"].loc[ARGDICT["star_solo_chem"], "barcode"]
        umi_len = ARGDICT["chemistry_file"].loc[ARGDICT["star_solo_chem"],"UMI"]
        
        if namey in starsolo_fastq_info:
        
            command = (f'STAR '
                       f'--runThreadN {ARGDICT["nthreads"]} '
                       f'--genomeDir {ARGDICT["star_solo_reference"]} '
                       f'--readFilesIn {starsolo_fastq_info[namey]["R2"]} {starsolo_fastq_info[namey]["R1"]} '
                       f'--soloType CB_UMI_Simple '
                       f'--soloCBwhitelist {ARGDICT["barcodes"]}/{barcode_file} '
                       f'--readFilesCommand zcat '
                       f'--soloCBlen {barcode_len} '
                       f'--soloUMIlen {umi_len} '
                       f'--outFileNamePrefix {staroutputdir} '
                       f'--soloFeatures {ARGDICT["star_solo_features"]}')
            
            if "star_solo_bam" in ARGDICT:
                command += (' --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM '
                            '--outSAMtype BAM SortedByCoordinate')
            
            print(command)
            subprocess.run(command, shell=True)
        
        else:
            print(f'Sample {namey} is not available! Please check the name and retry')
            
    print("ORGANIZING SOLO OUTPUT")
    organize_starsolo_output(ARGDICT, samples)
   
   
# for multiprocessing
def execute_sbp_command(command, 
                        workdir='/data/cr_count_organized_output/'):
    '''execute the command using subprocess'''
    try:
        output = subprocess.check_output(command,
                                         shell=True,
                                         cwd=workdir,
                                         stderr=subprocess.STDOUT)
        return output.decode()
    except subprocess.CalledProcessError as errata:
        return errata.output.decode()


def run_alignment_qc(ARGDICT, alignment_type):
    '''run alignment qc on cellranger count bams'''
    
    if alignment_type == 'cr_count':
    
        alignment_dir = '/data/cr_count_organized_output/cr_count_alignment_qc'
        
        if not os.path.isdir(alignment_dir):
            os.mkdir(alignment_dir)
        
        bams = sorted(glob.glob('/data/cr_count_organized_output/bam_files/*.bam'))
        ave_bam_load = sum([os.path.getsize(bam) / 2**30 for bam in bams])/len(bams)
        samples = [os.path.basename(bam).split('_possorted')[0] for bam in bams]
    
    if alignment_type == 'star':
        
        alignment_dir = '/data/STAR_organized_output/star_alignment_qc'
        
        if not os.path.isdir(alignment_dir):
            os.mkdir(alignment_dir)
        
        bams = sorted(glob.glob('/data/STAR_output/*.bam'))
        ave_bam_load = sum([os.path.getsize(bam) / 2**30 for bam in bams])/len(bams)
        
        # could in the future be multiprocessed
        for bam in bams:
            if not os.path.isfile(f'{bam}.bai'):
                pysam.index(bam)
        
        samples = [os.path.basename(bam).split('Aligned')[0] for bam in bams] 
    
    # some functions have high memory load
    # which can crash the program during multiprocessing
    # this part will adjust the number of threads
    total_bam_load = ave_bam_load * ARGDICT['nthreads']
    actual_bam_load = total_bam_load*0.15 # very approx
    mem = psutil.virtual_memory()
    available_memory_gb = mem.available / 2**30
    if actual_bam_load > available_memory_gb:
        revised_threads = int(available_memory_gb*(1/0.15) / ave_bam_load)
    else:
        revised_threads = ARGDICT['nthreads']

    for sample in samples:
        if not os.path.isdir(f'{alignment_dir}/{sample}'):
            os.mkdir(f'{alignment_dir}/{sample}')
            #os.mkdir(f'{alignment_dir}/{sample}/work')

    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f'performing GC stats at {current_time}')
    commands = [(f'read_GC.py '
                 f'-i {bam} '
                 f'-o {sample}',
                 f'{alignment_dir}/{sample}')
                 for bam, sample in zip(bams, samples)]
 
    with Pool(processes=ARGDICT["nthreads"]) as pool:
        pool.starmap(execute_sbp_command, commands)
    
    # convert to bam
    anno_dir = '/data/cr_count_reference/genes'
    if os.path.isfile(f'{anno_dir}/genes.gtf.gz'):
        gff_to_bed = (f'zcat {anno_dir}/genes.gtf.gz | bedparse gtf2bed \
                      > {anno_dir}/genes.bed')
        subprocess.run(gff_to_bed, shell=True)
        print('conversion of GTF to BED12 a success!')
    elif os.path.isfile(f'{anno_dir}/genes.gtf'):
        gff_to_bed = (f'bedparse gtf2bed {anno_dir}/genes.gtf \
                      > {anno_dir}/genes.bed')
        print('conversion of GTF to BED12 a success!')
        subprocess.run(gff_to_bed, shell=True)
    else:
        sys.exit('GTF not found; cannot convert to BED. exiting.')
    
    if "alignment_qc_genebody" in ARGDICT:
        
        now = datetime.now()
        current_time = now.strftime("%Y-%m-%d %H:%M:%S")
        print(f'performing gene body coverage at {current_time}')
        
        commands = [(f'geneBody_coverage.py -i {bam} '
                     f'-r {anno_dir}/genes.bed '
                     f'-o {sample}',
                     f'{alignment_dir}/{sample}')
                    for bam, sample in zip(bams, samples)]
    
        with Pool(processes=revised_threads) as pool:
            pool.starmap(execute_sbp_command, commands)
    
    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f'performing junction annotation at {current_time}')
    
    commands = [(f'junction_annotation.py -i {bam} '
                 f'-r {anno_dir}/genes.bed '
                 f'-o {sample}',
                 f'{alignment_dir}/{sample}')
                for bam, sample in zip(bams, samples)]

    with Pool(processes=ARGDICT["nthreads"]) as pool:
        pool.starmap(execute_sbp_command, commands)
        
    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f'performing junction saturation at {current_time}')
        
    commands = [(f'junction_saturation.py -i {bam} '
                 f'-r {anno_dir}/genes.bed '
                 f'-o {sample}',
                 f'{alignment_dir}/{sample}')
                for bam, sample in zip(bams, samples)]

    with Pool(processes=revised_threads) as pool:
        pool.starmap(execute_sbp_command, commands)
        
    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f'performing read distribution at {current_time}')
    
    print('performing read distribution')
    commands= [(f'read_distribution.py '
                f'-i {bam} '
                f'-r {anno_dir}/genes.bed '
                f'> {sample}.read_dist.txt',
                f'{alignment_dir}/{sample}')
                     for bam, sample in zip(bams, samples)]
                
    with Pool(processes=ARGDICT["nthreads"]) as pool:
        pool.starmap(execute_sbp_command, commands)
    
    print('alingment QC done')


def run_alignment_qc_master(ARGDICT):
    '''run alignment qc on cellranger count and star solo'''
    
    if not os.path.isdir('/data/cr_count_organized_output/bam_files'):
        print('CELLRANGER QC NOT PERFORMED BECAUSE BAMS NOT GENERATED')
    else:
        run_alignment_qc(ARGDICT, alignment_type = "cr_count")
    
    if not os.path.isdir('/data/STAR_output') or len(glob.glob('/data/STAR_output/*.bam'))==0:
        print('STAR QC NOT PERFORMED BECAUSE BAMS NOT GENERATED')
    else:
        run_alignment_qc(ARGDICT, alignment_type = "star")

        
def run_multi_qc(ARGDICT):
    '''this function runs multiqc'''
    
    multiqc_outdir = '/data/multiqc_report'
    
    if not os.path.isdir(multiqc_outdir):
        os.mkdir(multiqc_outdir)
    
    command = [rsdir for rsdir in ['/data/fastqc_output/',
                                    '/data/cr_count_organized_output/',
                                    '/data/STAR_organized_output/']
                if os.path.isdir(rsdir)]
    
    if len(command) > 0:
        command = ['multiqc'] + command
        command = ' '.join(command)
    
        subprocess.run(command, shell=True,
                       cwd=multiqc_outdir)