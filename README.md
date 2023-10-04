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

The pipeline inputs (and for that matter, outputs) are all contained in single folder, hereafter named ```workdir``` (but can be named whatever you'd like).
How to name the <run_id> folders is up to you. We recommend using something recognizable like the flow cell number.

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

Place the reference genome in the ```cr_count_reference``` folder, and your BCLs into the appropriate folders:

```
data
├── input_bcl
│   ├── <run_id>
│   └── <run_id>
└── cr_count_reference
```

### I have FASTQs from two sequencing runs and I want to get count matrices using STARsolo.

Place the reference genome in the ```star_solo_reference``` folder, and your BCLs into the appropriate folders:

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

### I have BCLs from a new sequencing run, FASTQs from a previous run, and I would like to get count matrices using Cellranger. However, I need to modify the human reference genome with a custom GFP gene used in my experiment.

Place the ```genome.fa``` and ```genes.gtf``` in the ```cr_count_reference_template``` folder (Note: this is will also work if you place them into the ```star_solo_reference_template``` folder for STARsolo), and your BCLs/FASTQs into the appropriate folders:

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
└── cr_count_reference_template
    ├── genome.fa
    ├── genes.gtf
    └── add.fa
```

The ```add.fa``` file should be a simple FASTA file, for example (and as per Cell Ranger):

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

## Perform QC

<details>
<br>

The pipeline is equipped to run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on FASTQ files by using the flag ```--fastqc```

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
    journal = {Journal},
    year = {2023},
    url = {URL},
    doi = {DOI},
  }
```
</details>



