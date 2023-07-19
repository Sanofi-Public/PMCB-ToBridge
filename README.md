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


## Quick Start
Want to quickly convert your binary base call (BCL) to FASTQs, run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and then align them to your reference genome with [cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) and get neatly organized output that can be directly analyzed with CellBridge?

This simple command (given that you have the Docker image built or pulled) will do the trick:

```docker run -it -v /root/workdir/:/data:z tobridge tobridge --project my2bridge_run --bcl_convert --cr_count --fastqc```

**Important:** this will run in your working directory (in this case */root/workdir*). Your BCL files and the corresponding SampleSheet.csv in the proper format must be (within the working directory) in the directory called *input_bcl* and your reference genome must be in the directory called *cr_count_reference*. For cellranger count, you can start by downloading the pre-built reference genomes from [10x](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) as tar archives.

## Cellranger mkfastq and STARsolo
Additionally, ToBridge supports BCL conversion using [cellranger mkfastq](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq), which is slower (but comes with 10x support and the ability to use 10x indices in SampleSheet.csv) and alignment using [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md), which is considerably faster than cellranger count, but does not come with 10x support and requires more information to run. To run mkfastq, use the ```--cr_mkfastq``` flag. To run STARsolo, you should know the correct chemistry, for example,  ```--star_solo --star_solo_chem SC3Pv3``` for Single Cell 3' v3. Reference genome will need to be in the *star_solo_reference* directory and should only include the STAR indices. An option to build your own reference genome will be added to a later version of ToBridge. For now, please contact us for a human reference if you do not have one.

## Additional Options
You can provide additional FASTQ files in the *input_fastq* folder. Since the same library is often sequenced multiple times, we have included the ability to align multiple rounds of sequencing. These FASTQ files can provided in the folder directly or in subfolders within (for example, you can place two sets of FASTQs in folders *previous_fastq1* and *previous_fastq2*, and place them in the *input_fastq* folder. Additionally, if you would only like to align a select number of samples for the time being, you can provide an argument e.g, ```--samp S1,S5```

## Contact
For help and questions please contact the [cellbridge's maintenance team](mailto:nima.nouri@sanofi.com).
