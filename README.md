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

Place the reference genome in the ```cr_count_reference``` folder, and your BCLs into the appropriate folders:

```
workdir
├── input_bcl
│   ├── <run_id>
│   └── <run_id>
└── cr_count_reference
```

### I have FASTQs from two sequencing runs and I want to get count matrices using STARsolo.

Place the reference genome in the ```star_solo_reference``` folder, and your BCLs into the appropriate folders:

```
workdir
├── input_fastq
│   ├── <run_id>
│   │    ├── S1_S1_L001_R1_001.fastq.gz
│   │    ├── S1_S1_L001_R1_001.fastq.gz
│   │    ├── S2_S2_L001_R1_001.fastq.gz
│   │    └── S2_S2_L001_R1_001.fastq.gz
│   └── <run_id>
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S2_S2_L001_R1_001.fastq.gz
│        └── S2_S2_L001_R1_001.fastq.gz
└── star_solo_reference
```

### I have BCLs from a new sequencing run, FASTQs from a previous run, and I would like to get count matrices using Cellranger.
### However, I need to modify the human reference genome with a custom GFP gene used in my experiment.

Place the ```genome.fa``` and ```genes.gtf``` in the ```cr_count_reference_template``` folder (Note: this is will also work if you place them into the ```star_solo_reference_template``` folder for STARsolo), and your BCLs/FASTQs into the appropriate folders:

```
workdir
├── input_bcl
│   └── <run_id>
├── input_fastq
│   └── <run_id>
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S1_S1_L001_R1_001.fastq.gz
│        ├── S2_S2_L001_R1_001.fastq.gz
│        └── S2_S2_L001_R1_001.fastq.gz
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

</details>

---

