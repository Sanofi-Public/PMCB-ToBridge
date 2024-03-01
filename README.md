[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

<p align="center" width="100%">
<img width="40%" src="./cellbridge_logo.png"> 
</p>

# News

<b>Feb 8, 2024</b>: major updates to QC in <b>ToBridge v1.1.0</b>!
We have incorporated alignment QC as well as comprehensive QC reporting with [MultiQC](https://multiqc.info/).

We are constantly looking for ways to improve the pipeline. Currently on our radar:

--> [velocyto](https://velocyto.org/velocyto.py/tutorial/cli.html#run10x-run-on-10x-chromium-samples) on 10X samples

--> simplifying how feature barcoding is setup

--> an informative log file

# CellBridge

**CellBridge** is an automated and versatile workflow meticulously designed to
simplify and expedite the standard procedures entailed in scRNA-seq analysis,
eliminating the need for specialized bioinformatics expertise. CellBridge
harnesses cutting-edge computational methods, integrating an array of advanced
functionalities. It encompasses various crucial steps in scRNA-seq analysis,
starting from the initial conversion of raw unaligned sequencing reads into the
FASTQ format, followed by read alignment, gene expression quantification,
quality control, doublet removal, normalization, batch correction,
dimensionality reduction, clustering, identification of cell markers, and
accurate cell type annotation. CellBridge provides convenient parameterization
of the workflow, while its Docker-based framework ensures reproducibility of
results across diverse computing environments. 

The CellBridge ecosystem comprises two main executables: 1) `tobridge` for pre-processing and 2) `cellbridge` for processing. Please see [CellBridge](https://github.com/Sanofi-Public/PMCB-CellBridge) Github page for the processing step of the data.

<p align="center" width="100%">
<img width="85%" src="./pipeline_schematic.png"> 
</p>

---

## Workflow Inputs

<details>

The pipeline inputs (and for that matter, outputs) are all contained in single folder, hereafter named ```workdir``` (but can be named whatever you'd like).
How to name the <run_id> folders is up to you. We recommend using something recognizable like the flow cell number. 

Each BCL folder should contain a ```SampleSheet.csv```. Please refer to [10X](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq#simple_csv) if you are in doubt on how to create one.


```
data
├── input_bcl
│   ├── <run_id>
│   └── <run_id>
├── input_fastq
│   ├── <run_id>
│   └── <run_id>
├── cr_count_reference
├── star_solo_reference
├── cr_count_reference_template
│   ├── genome.fa
│   └── genes.gtf
├── star_solo_reference_template
│   ├── genome.fa
│   └── genes.gtf
├── library_files 
└── feature_ref.csv
```

### I have BCLs form two sequencing runs and I want to get count matrices using Cellranger count.

Place the Cell Ranger reference genome in the ```cr_count_reference``` folder, and your BCLs into the appropriate folders:

```
data
├── input_bcl
│   ├── <run_id>
│   └── <run_id>
└── cr_count_reference
```

### I have FASTQs from two sequencing runs and I want to get count matrices using STARsolo.

Place the STARsolo reference genome in the ```star_solo_reference``` folder, and your BCLs into the appropriate folders:

```
data
├── input_fastq
│   ├── <run_id>
│   │    ├── S1_S1_L001_R1_001.fastq.gz
│   │    ├── S1_S1_L001_R2_001.fastq.gz
│   │    ├── S2_S2_L001_R1_001.fastq.gz
│   │    └── S2_S2_L001_R2_001.fastq.gz
│   └── <run_id>
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S1_S1_L001_R2_001.fastq.gz
│        ├── S2_S2_L001_R1_001.fastq.gz
│        └── S2_S2_L001_R2_001.fastq.gz
└── star_solo_reference
```

### I have BCLs from a new sequencing run, FASTQs from a previous run, and I would like to get count matrices using STARsolo. However, I need to modify the human reference genome with a custom GFP gene used in my experiment.

Place the ```genome.fa``` and ```genes.gtf``` files in the ```star_solo_reference_template``` folder (Note: this is will also work if you place them into the ```cr_count_reference_template``` folder for Cell Ranger), and your BCLs/FASTQs into the appropriate folders. If you do not have readily available genome.fa and genes.gtf, you can use ones from a publicly available [10X repository](https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads). Note that you can create a reference with ```genome.fa``` and ```genes.gtf``` alone; the ```add.fa``` is optional and comes in handy if you want to add sequences to both your genome and the GTF annotation.

```
data
├── input_bcl
│   └── <run_id>
├── input_fastq
│   └── <run_id>
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S1_S1_L001_R2_001.fastq.gz
│        ├── S2_S2_L001_R1_001.fastq.gz
│        └── S2_S2_L001_R2_001.fastq.gz
└── star_solo_reference_template
    ├── genome.fa
    ├── genes.gtf
    └── add.fa
```

The ```add.fa``` file should be a simple FASTA file, for example (and as per an example provided by 10X):

```
>GFP
TACACACGAATAAAAGATAACAAAGATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTT
GTTGAATTAGATGGCGATGTTAATGGGCAAAAATTCTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACAT
ACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGGAAGCTACCTGTTCCATGGCCAACACTTGTCAC
TACTTTCTCTTATGGTGTTCAATGCTTTTCAAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAG
AGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTACAAAGATGACGGGAACTACAAGACAC
GTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGA
AGATGGAAACATTCTTGGACACAAAATGGAATACAACTATAACTCACATAATGTATACATCATGGCAGAC
AAACCAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTAAAGATGGAAGCGTTCAATTAG
CAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTC
CACACAATCTGCCCTTTCCAAAGATCCCAACGAAAAGAGAGATCACATGATCCTTCTTGAGTTTGTAACA
GCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAAATGTCCAGACTTCCAATTGACACTAAAG
TGTCCGAACAATTACTAAATTCTCAGGGTTCCTGGTTAAATTCAGGCTGAGACTTTATTTATATATTTAT
AGATTCATTAAAATTTTATGAATAATTTATTGATGTTATTAATAGGGGCTATTTTCTTATTAAATAGGCT
ACTGGAGTGTAT
```

### I have FASTQs of CITE-seq data that I would like to process using Cell Ranger.

Place the reference genome in the ```cr_count_reference``` folder, and your FASTQs into the appropriate folder. 
Include the feature reference sequences in  ```feature_ref.csv``` as per [Cell Ranger instructions](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis), and place your sample descriptions in the ```library_files``` folder (more on this below):

```
data
├── input_fastq
│   └── <run_id>
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S1_S1_L001_R2_001.fastq.gz
│        ├── S2_S2_L001_R1_001.fastq.gz
│        ├── S2_S2_L001_R2_001.fastq.gz
│        ├── S3_S3_L001_R1_001.fastq.gz
│        ├── S3_S3_L001_R2_001.fastq.gz
│        ├── S4_S4_L001_R1_001.fastq.gz
│        └── S4_S4_L001_R2_001.fastq.gz
├── cr_count_reference
├── library_files 
└── feature_ref.csv
```

Format your ```library_files``` as per [Cell Ranger instructions](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis).
For example, if S1 and S2 come from your antibody capture libraries and S3 and S4 come from the corresponding gene expression libraries, you would make two files as follows:

```S3_S1.csv```:

```
fastqs,sample,library_type
/data/input_fastq/<run_id>,S1,Antibody Capture
/data/input_fastq/<run_id>,S3,Gene Expression
```

```S4_S2.csv```:

```
fastqs,sample,library_type
/data/input_fastq/<run_id>,S2,Antibody Capture
/data/input_fastq/<run_id>,S4,Gene Expression
```

and place them in the ```library_files``` directory.

</details>

---

## Docker Images

<details>
<br>

The pre-built images are available in the `pmcbscb` (Precision Medicine and
Computational Biology – Single Cell Biology) Docker Hub repository. They can be
seamlessly pulled by:

```
docker pull pmcbscb/tobridge
```
```
docker pull pmcbscb/cellbridge
```

Note: Please see [CellBridge](https://github.com/Sanofi-Public/PMCB-CellBridge) for the main
processing steps.

</details>

---

## Flag Options

<details>
<br>

The extensive documentation for flag options is embedded within the workflows.
For a review of the flags, please execute:

```
docker run pmcbscb/tobridge tobridge --help
```
```
docker run pmcbscb/cellbridge cellbridge --help
```

For detailed information about the available flag options, refer
to our up-to-date HTML manual:
[tobridge-flags](http://htmlpreview.github.io/?https://github.com/Sanofi-Public/PMCB-ToBridge/blob/master/tobridge_flags.html)

Note: Please see [CellBridge](https://github.com/Sanofi-Public/PMCB-CellBridge) for the main
processing steps.

</details>

---

## Demo Workflow

<details>
<br>

#### Get fastq demo files

Users can download FASTQ files from one of the publicly-available data sets on
the 10x Genomics support site. This example uses the 1,000 PBMC data set from
human peripheral blood mononuclear cells (PBMC), consisting of lymphocytes (T
cells, B cell, and NK kills) and monocytes. Please copy and paste the following
instructions into the terminal:

``` 
mkdir sandbox && cd sandbox && \
mkdir -p input_fastq/run_1 && \
wget -P input_fastq/run_1 https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar && \
tar -xvf input_fastq/run_1/pbmc_1k_v3_fastqs.tar -C input_fastq/run_1 --strip-components=1
```

#### Get the reference transcriptome and metadata

The following command lines set up the required data structure to run the workflow:

``` 
mkdir cr_count_reference && \
wget -P cr_count_reference https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz && \
tar -zxvf cr_count_reference/refdata-gex-GRCh38-2020-A.tar.gz -C cr_count_reference --strip-components=1
```

#### Execute workflows

Assuming the images have already been pulled (see 'Docker Images' above), run
the pre-processing pipeline which perform `fastqc` and `cellranger count`:

``` 
docker run -v ${PWD}:/data:z pmcbscb/tobridge:latest tobridge \
                                           --fastqc \
                                           --cr_count \
                                           --cr_count_bam \
                                           --alignment_qc
```

Move to `cellbridge` directory and download the `metadata`:

```
cd cr_count_organized_output/cellbridge && \
wget https://raw.githubusercontent.com/Sanofi-Public/PMCB-CellBridge/master/demo/metadata.csv 
```

Run the processing pipeline with the basic flags:

``` 
docker run -v ${PWD}:/data:z pmcbscb/cellbridge:latest cellbridge \
                                           --project project-demo \
                                           --species hs \
                                           --tissue pbmc \
                                           --metadata sample_based
```

Note: sharing files between the host operating system and the container requires
you to bind a directory on the host to the container mount points using the `-v`
argument. There is one available mount points defined in the container named
`data`. In the example above the current directory `${PWD}` was used and not an
absolute notation. If you intended to pass a host directory, use absolute path.

Note: for details about the processing step, please visit the
[CellBridge](https://github.com/Sanofi-Public/PMCB-CellBridge) GitHub page.

</details>

---

## Workflow Outputs

<details>
<br>

The <b>pre-processing</b> part of the pipeline (i.e. ToBridge) has the following output structure:

```
data
├── STAR_organized_output
│   ├── cellbridge
│   ├── metrics.csv
│   └── star_alignment_qc
├── STAR_output
├── cr_count_organized_output
│   ├── cellbridge
│   ├── cr_count_alignment_qc
│   ├── loupe_files
│   ├── metrics.csv
│   └── web_summaries
├── cr_count_output
├── fastqc_output
└── multiqc_report
```

While some of these are self-explanatory, others call for additional clarification.

```cellbridge``` directories have the folder structure ready to be plugged into the main portion of the pipeline.

In ```cr_count_organized_output```, Loupe files and web summaries are grouped together for all the samples, and ```metrics.csv``` has the metrics for all the samples in the same file.
Ditto for ```STAR_organized_output``` with respect to ```metrics.csv```.

Raw [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) outputs and [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) outputs are found in ```STAR_output``` and ```cr_count_output```, respectively.

[RSeQC](https://rseqc.sourceforge.net/) alignment QC for STARsolo and Cell Ranger is outputted in subdirectories ```star_alignment_qc``` and ```cr_count_alignment_qc```, respectively.

The combined QC output by [MultiQC](https://multiqc.info/) is in the ```multiqc_report``` folder.

Note: for outputs of the processing step, please visit
[CellBridge](https://github.com/Sanofi-Public/PMCB-CellBridge) Github page."

</details>

---

## Perform QC

<details>
<br>

The pipeline is equipped to run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on FASTQ files by using the flag ```--fastqc```

It also provides [RSeQC](https://rseqc.sourceforge.net/)-based alignment QC on BAM files.
(please note that to get BAMs to output with Cell Ranger, you'll need to use the ```--cr_count_bam``` flag).
Specifically, it will determine GC count, [junction annotation](https://rseqc.sourceforge.net/#junction-annotation-py) and saturation, and [read distribution](https://rseqc.sourceforge.net/#read-distribution-py).
This QC will run if you use the the ```--alignment_qc``` flag. 

If you would like information on [gene body coverage](https://rseqc.sourceforge.net/#genebody-coverage-py) as well, please use the flag ```--alignment_qc_genebody``` instead.
Please note that gene body coverage calculations are time-consuming.

The QC report combining FastQC, RSeQC, and the alignment QC reports automatically generated by Cell Ranger is produced by [MultiQC](https://multiqc.info/).

</details>

---

## Contact

<details>
<br>

For help and questions please contact the [cellbridge's maintenance team](mailto:nima.nouri@sanofi.com).

</details>

---

## Citing CellBridge

<details>
<br>

If you use CellBridge please cite our paper: 

```
  @Article{,
    author = {Nima Nouri and Andre H. Kurlovs, et al.},
    title = {Scaling up Single-Cell RNA-seq Data Analysis with CellBridge Workflow},
    journal = {Bioinformatics},
    year = {2023},
    url = {https://academic.oup.com/bioinformatics/article/39/12/btad760/7479685},
    doi = {10.1093/bioinformatics/btad760},
  }
```
</details>
