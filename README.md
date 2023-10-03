[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

<p align="center" width="100%">
<img width="40%" src="./cellbridge_logo.png"> 
</p>

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

See [CellBridge](https://github.com/Sanofi-Public/PMCB-CellBridge) for the
processing of the data.

<p align="center" width="100%">
<img width="85%" src="./pipeline_schematic.png"> 
</p>

---

## Workflow Inputs

<details>
<br>

The pipeline inputs (and for that matter, outputs) are all contained in single folder, hereafter named workdir (but can be named whatever you would like).

```
workdir
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

Which folders you include depends on your specific case. Some examples (from simpler to more complex) are included below:

### I have BCLs form two sequencing runs and I want to get count matrices using Cellranger count.

Place the reference genome in the cr_count_reference_template folder, and your BCLs into the approproate folders:
```
workdir
├── input_bcl
│   ├── <run_id>
│   └── <run_id>
└── cr_count_reference
```

### I have FASTQs from two sequencing runs and I want to get count matrices using STARsolo.

```
workdir
├── input_fastq
│   ├── <run_id>
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S2_S2_L001_R1_001.fastq.gz
│        └── S2_S2_L001_R1_001.fastq.gz
│   └── <run_id>
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S2_S2_L001_R1_001.fastq.gz
│        └── S2_S2_L001_R1_001.fastq.gz
└── star_solo_reference
```

### I have BCLs from a new sequencing run, FASTQs from a previous run, and I would like to get count matrices using Cellranger.
### However, I need to modify the human reference genome with a custom GFP gene used in my experiment.

```
workdir
├── input_bcl
│   └── <run_id>
├── input_fastq
│   └── <run_id>
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S2_S2_L001_R1_001.fastq.gz
│        ├── S2_S2_L001_R1_001.fastq.gz
├── cr_count_reference_template
│   ├── genome.fa
│   └── genes.gtf


```




</details>

---


## Quick Start
Want to quickly convert your binary base call (BCL) to FASTQs, run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and then align them to your reference genome with [cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) and get neatly organized output that can be directly analyzed with CellBridge?

This simple command (given that you have the Docker image built or pulled) will do the trick:

```docker run -it -v /root/workdir/:/data:z <docker_image_name> tobridge --bcl_convert --cr_count```

**Important, Part 1:** a public pre-built Docker image is provided! Please see the Docker section below.

**Important, Part 2:** this will run in your working directory (in this case */root/workdir*). Your BCL files and the corresponding SampleSheet.csv must be (within the working directory) in the directory called *input_bcl* and your reference genome must be in the directory called *cr_count_reference*. If you have a 10X-formatted sample sheet, please also include the flag ```--bcl_convert_sheet_conv```. For cellranger count, you can start by downloading the pre-built reference genomes from [10x](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) as tar archives. The FASTQ files will end up in the *input_fastq* folder (named as such because they will be used as inputs to the alignment.

The raw output will be in the *cr_count_output* directory while a more readable and organized output will be placed in the *cr_count_organized_output* directory.

## Cell Ranger's mkfastq and STARsolo
Additionally, ToBridge supports BCL conversion using [cellranger mkfastq](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq), which is slower (but comes with 10x support) and alignment using [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md), which is considerably faster than cellranger count, but does not come with 10x support and requires more information to run. To run mkfastq, use the ```--cr_mkfastq``` flag. To run STARsolo, you should know the correct chemistry, for example,  ```--star_solo --star_solo_chem SC3Pv3``` for Single Cell 3' v3. 

Reference genome will need to be in the *star_solo_reference* directory and should only include the STAR indices. For your convenience, we provide an option to build your own reference genome with the flag ```--star_solo_genome_generate``` based on user-provided *genome.fa* and *genes.gtf*, both of which need to be placed in the *star_solo_reference_template* directory. If in doubt, you can obtain these two files from the Cell Ranger reference build!

## Custom reference genomes
You can build your own reference genome by putting your custom *genome.fa* and *genes.gtf* in the *star_solo_reference_template* for STAR Solo and in the *cr_count_reference_template* directory for Cell Ranger. Since it's tedious to edit the GTF, we have introduced a way to do it by simply adding *add.fa* to the folder(s). The *add.fa* file should be a standard FASTA file with the sequences you would like to add to your reference genome in the + orientation (i.e. please reverse-complement beforehand them if needed). The pipeline will automatically update both the *genome.fa* and *genes.gtf* based on *add.fa* prior to generating a new reference genome.

## Feature barcode analysis with Cell Ranger
We have incorporated the ability to analyze feature barcodes with [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis). Simply place your *feature_ref.csv* in the working directory, your library file CSVs in the *library_files* directory, and provide the flag ```--cr_count_feature```.

## Additional Options
You can provide additional FASTQ files in the *input_fastq* folder. Since the same library is often sequenced multiple times, we have included the ability to align multiple rounds of sequencing. For example, you can place two sets of previously convereted FASTQ files in folders *previous_fastq1* and *previous_fastq2*, and place them in the *input_fastq* folder. Moreover, you can combine BCL conversion of a new sequencing run with previously converted FASTQ files.

If you would only like to align a select number of samples for the time being, you can provide an argument e.g, ```--samp S1,S5```

## Docker
If you would like to build the package from scratch, please note that you will have to download some of the files yourself since 10X and Illumina require a login and provide a different key each time. Please refer to the Dockerfile for more details. A public Docker image of this pipeline is temporarily (i.e. until it finds a more permanent home) available from one of the developers' Docker repos by running:

```docker pull kurlovs/cellbridge:tobridge.v0.2.4```

To use this public Docker image after pulling it, specify its name in the command. For example:

```docker run -it -v /root/workdir/:/data:z kurlovs/cellbridge:tobridge.v0.2.4 tobridge --bcl_convert --bcl_convert_sheet_conv --cr_count```

## Contact
For help and questions please contact the [cellbridge's maintenance team](mailto:nima.nouri@sanofi.com).
