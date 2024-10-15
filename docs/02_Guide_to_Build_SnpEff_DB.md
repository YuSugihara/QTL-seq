# Building a Custom SnpEff Database

## Table of Contents
1. [Introduction](#1-introduction)  
2. [Downloading Reference and Annotation Data](#2-downloading-reference-and-annotation-data)  
3. [Setting Paths and Navigating to the snpEff Directory](#3-setting-paths-and-navigating-to-the-snpeff-directory)  
4. [Editing snpEff.config](#4-editing-snpeffconfig)  
5. [Building the SnpEff Database](#5-building-the-snpeff-database)  
6. [Handling Errors During Build](#6-handling-errors-during-build)  
7. [Verifying Mutations with IGV](#7-verifying-mutations-with-igv)

---

## 1. Introduction

In this tutorial, we will guide you through the process of building a custom SnpEff database using **NCBI RefSeq Annotations** for **rice (Oryza sativa)**. The steps outlined here are meant for users who need to annotate genomic variants using their own reference genome or custom annotations that may not be available in the pre-built SnpEff databases.

Before proceeding, it is **highly recommended** to carefully read the [SnpEff documentation](https://pcingola.github.io/SnpEff/snpeff/introduction/). This tutorial assumes that you have already installed SnpEff using Conda (`conda install -c bioconda snpeff`). If your reference genome is already available in a pre-built SnpEff database, you can avoid building a custom database. For QTL-seq users, please refer to the [QTL-seq manual](https://github.com/YuSugihara/QTL-seq) for details on how to use the `-e` option to select pre-built databases.

## 2. Downloading Reference and Annotation Data

For this tutorial, we will use the **NCBI RefSeq Annotations** for **rice (Oryza sativa)**. The dataset can be accessed at the following URL:

[NCBI RefSeq Annotations for Oryza sativa](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_034140825.1/)

We will use the following links to download the reference FASTA and GFF files using `wget`, then decompress them using `gunzip`:

- Reference FASTA: [GCF_034140825.1_ASM3414082v1_genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.fna.gz)
- GFF file: [GCF_034140825.1_ASM3414082v1_genomic.gff.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.gff.gz)

To download and decompress these files, run the following commands:

```bash
# Download reference FASTA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.fna.gz

# Download GFF file
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.gff.gz

# Decompress the files
gunzip GCF_034140825.1_ASM3414082v1_genomic.fna.gz
gunzip GCF_034140825.1_ASM3414082v1_genomic.gff.gz
```

Next, define the paths to these files as variables:

```bash
PATH_TO_FASTA=$(pwd)/GCF_034140825.1_ASM3414082v1_genomic.fna
PATH_TO_GFF=$(pwd)/GCF_034140825.1_ASM3414082v1_genomic.gff
```

## 3. Setting Paths and Navigating to the snpEff Directory

Once the reference and annotation files are downloaded and paths are set, you need to navigate to the directory where `snpEff` is installed. If you installed `snpEff` using Conda, it is typically located at:

```
${HOME}/anaconda/envs/qtlseq/share/snpeff-5.2-1
```

However, the exact path may vary depending on whether you are using Conda or Mamba, and it may also differ based on the version of `snpEff` you have installed.

After navigating to this directory, you should see the following directories and files:

- `scripts`
- `snpEff`
- `snpEff.config`
- `snpEff.jar`

Now, create a directory for your custom genome data:

```bash
mkdir -p data/GCF_034140825.1
```

Next, create symbolic links for the reference FASTA and GFF files in the newly created directory:

```bash
ln -s ${PATH_TO_FASTA} data/GCF_034140825.1/sequences.fa
ln -s ${PATH_TO_GFF} data/GCF_034140825.1/genes.gff
```

## 4. Editing snpEff.config

After setting up the files, the next step is to edit the `snpEff.config` file to register your custom genome. Open `snpEff.config` in a text editor and add the following line:

```
GCF_034140825.1.genome : GCF_034140825.1
```

Make sure that the second `GCF_034140825.1` matches the directory name you created in the `data/` folder. This line tells `snpEff` where to find the genome data for your custom annotation.

Save the file after adding the line.

## 5. Building the SnpEff Database

Once the `snpEff.config` file is updated, you can start building the SnpEff database using the following command:

```bash
snpEff build -gff3 -noCheckCds -noCheckProtein GCF_034140825.1
```

This command builds the database using the reference genome and annotation files located in `data/GCF_034140825.1/`.

## 6. Handling Errors During Build

If the GFF file is correctly formatted and compatible with `snpEff`, the build process will complete successfully, and you will be able to annotate variants using the `-e` option with the specified genome.

However, in some cases, you might encounter errors such as:

```
WARNING_CHROMOSOME_CIRCULAR: Chromosome 'NC_001320.1'
WARNING_GENE_NOT_FOUND: Gene 'null' (NC_001751.1:1732-2019)
```

These indicate that the build process has failed due to issues with annotations for certain contigs (e.g., organelle sequences) in the GFF file.

### Removing Problematic Contigs from the GFF

To fix this issue, you can remove the problematic contig annotations from the GFF file before building the database. Below is an example of how to do this using `grep` to filter out both `NC_001320.1` and `NC_001751.1`:

```bash
# Create a new GFF file excluding the problematic contigs
grep -v -e "NC_001320.1" -e "NC_001751.1" ${PATH_TO_GFF} > filtered_genes.gff

# Update the symbolic link to point to the new filtered GFF file
ln -sf $(pwd)/filtered_genes.gff data/GCF_034140825.1/genes.gff
```

This command will exclude any lines in the GFF file that contain `NC_001320.1` or `NC_001751.1`, the contigs causing the issue. After creating the new `filtered_genes.gff` file, update the symbolic link for the GFF file in the `snpEff` data directory.

### Re-running the Build Command

Once the GFF file has been filtered, you can re-run the `snpEff build` command:

```bash
snpEff build -gff3 -noCheckCds -noCheckProtein GCF_034140825.1
```

If the problematic contigs were the only cause of the error, the build should now complete successfully.

## 7. Verifying Mutations with IGV

Even if your SnpEff database is built without errors, it is important to manually verify the annotations. GFF format inconsistencies may still cause SnpEff to annotate variants incorrectly, even when no errors are reported during the build process.

To check the accuracy of the annotations, you can use a genome visualization tool like IGV (Integrative Genomics Viewer). IGV allows you to visually inspect the variant locations and assess their impact on gene function.

You can download IGV and learn more about it from the following link:

[IGV](https://igv.org)
