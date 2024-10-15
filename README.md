# QTL-seq User Guide
#### version 2.2.7

## Table of contents
- [What is QTL-seq?](#what-is-qtl-seq)
- [Installation](#installation)
  + [Dependencies](#dependencies)
  + [Installation via bioconda](#installation-via-bioconda)
  + [Manual installation](#manual-installation)
- [Usage](#usage)
  + [Example 1 : run QTL-seq from FASTQ without trimming](#example-1--run-qtl-seq-from-fastq-without-trimming)
  + [Example 2 : run QTL-seq from FASTQ with trimming](#example-2--run-qtl-seq-from-fastq-with-trimming)
  + [Example 3 : run QTL-seq from BAM](#example-3--run-qtl-seq-from-bam)
  + [Example 4 : run QTL-seq from multiple FASTQs and BAMs](#example-4--run-qtl-seq-from-multiple-fastqs-and-bams)
  + [Example 5 : run QTL-plot from VCF](#example-5--run-qtl-plot-from-vcf)
- [Outputs](#outputs)
- [About multiple testing correction](#about-multiple-testing-correction)
- [Build and use your own database for snpEff](#build-and-use-your-own-database-for-snpeff)

## What is QTL-seq?
<img src="https://github.com/YuSugihara/QTL-seq/blob/master/images/1_logo.png" width=200>

Bulked segregant analysis, as implemented in  QTL-seq (Takagi et al., 2013), is a powerful and efficient method for identifying agronomically important loci in crop plants. QTL-seq has been adapted from MutMap to identify quantitative trait loci. It utilizes sequences pooled from two segregating progeny populations with extreme opposite traits (e.g. resistant vs susceptible) and a single whole-genome resequencing of either of the parental cultivars.

#### Citation
- Yu Sugihara, Lester Young, Hiroki Yaegashi, Satoshi Natsume, Daniel J. Shea, Hiroki Takagi, Helen Booker, Hideki Innan, Ryohei Terauchi, Akira Abe (2022). [High performance pipeline for MutMap and QTL-seq](https://doi.org/10.7717/peerj.13170). PeerJ, 10:e13170.

- Hiroki Takagi, Akira Abe, Kentaro Yoshida, Shunichi Kosugi, Satoshi Natsume, Chikako Mitsuoka, Aiko Uemura, Hiroe Utsushi, Muluneh Tamiru, Shohei Takuno, Hideki Innan, Liliana M. Cano, Sophien Kamoun, Ryohei Terauchi (2013).  [QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations](https://doi.org/10.1111/tpj.12105). Plant journal 74:174-183.


## Installation
### Dependencies
#### Software
- [BWA](http://bio-bwa.sourceforge.net/)
- [SAMtools](http://samtools.sourceforge.net/)
- [BCFtools](http://samtools.github.io/bcftools/)
- [SnpEff](http://snpeff.sourceforge.net/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

#### Python (>=3.5) libraries
- matplotlib
- numpy
- pandas
- seaborn (optional)

### Installation via bioconda
You can install QTL-seq via [bioconda](https://bioconda.github.io/index.html).
```
conda create -c bioconda -n qtlseq qtlseq
conda activate qtlseq
```

### Manual installation
If you encounter an error during installation, you can install QTL-seq manually.
```
git clone https://github.com/YuSugihara/QTL-seq.git
cd QTL-seq
pip install -e .
```
Then you need to install other dependencies yourself. We highly recommend installing SnpEff and Trimmomatic via bioconda.
```
conda install -c bioconda snpeff
conda install -c bioconda trimmomatic
```
After installation, please check whether SnpEff and Trimmomatic are working using the commands below.
```
snpEff --help
trimmomatic --help
```

## Usage

If your reference genome contains more than 50 contigs, only the 50 longest contigs will be plotted.

```
qtlseq -h

usage: qtlseq -r <FASTA> -p <BAM|FASTQ> -b1 <BAM|FASTQ>
              -b2 <BAM|FASTQ> -n1 <INT> -n2 <INT> -o <OUT_DIR>
              [-F <INT>] [-T] [-e <DATABASE>]

QTL-seq version 2.2.7

options:
  -h, --help         show this help message and exit
  -r , --ref         Reference FASTA file.
  -p , --parent      FASTQ or BAM file of the parent. If specifying
                     FASTQ, separate paired-end files with a comma,
                     e.g., -p fastq1,fastq2. This option can be
                     used multiple times.
  -b1 , --bulk1      FASTQ or BAM file of bulk 1. If specifying
                     FASTQ, separate paired-end files with a comma,
                     e.g., -b1 fastq1,fastq2. This option can be
                     used multiple times.
  -b2 , --bulk2      FASTQ or BAM file of bulk 2. If specifying
                     FASTQ, separate paired-end files with a comma,
                     e.g., -b2 fastq1,fastq2. This option can be
                     used multiple times.
  -n1 , --N-bulk1    Number of individuals in bulk 1.
  -n2 , --N-bulk2    Number of individuals in bulk 2.
  -o , --out         Output directory. The specified directory must not
                     already exist.
  -F , --filial      Filial generation. This parameter must be greater
                     than 1. [2]
  -t , --threads     Number of threads. If a value less than 1 is specified,
                     QTL-seq will use the maximum available threads. [2]
  -w , --window      Window size in kilobases (kb). [2000]
  -s , --step        Step size in kilobases (kb). [100]
  -D , --max-depth   Maximum depth of variants to be used. [250]
  -d , --min-depth   Minimum depth of variants to be used. [8]
  -N , --N-rep       Number of replicates for simulations to generate
                     null distribution. [5000]
  -T, --trim         Trim FASTQ files using Trimmomatic.
  -a , --adapter     FASTA file containing adapter sequences. This option
                     is used when "-T" is specified for trimming.
  --trim-params      Parameters for Trimmomatic. Input parameters
                     must be comma-separated in the following order:
                     Phred score, ILLUMINACLIP, LEADING, TRAILING,
                      SLIDINGWINDOW, MINLEN. To remove Illumina adapters,
                     specify the adapter FASTA file with "--adapter".
                     If not specified, adapter trimming will be skipped.
                     [33,<ADAPTER_FASTA>:2:30:10,20,20,4:15,75]
  -e , --snpEff      Predict causal variants using SnpEff. Check
                     available databases in SnpEff.
  --line-colors      Colors for threshold lines in plots. Specify a
                     comma-separated list in the order of SNP-index,
                     p95, and p99. ["#FE5F55,#6FD08C,#E3B505"]
  --dot-colors       Colors for dots in plots. Specify a
                     comma-separated list in the order of bulk1,
                     bulk2, and delta. ["#74D3AE,#FFBE0B,#40476D"]
  --mem              Maximum memory per thread when sorting BAM files;
                     suffixes K/M/G are recognized. [1G]
  -q , --min-MQ      Minimum mapping quality for mpileup. [40]
  -Q , --min-BQ      Minimum base quality for mpileup. [18]
  -C , --adjust-MQ   Adjust the mapping quality for mpileup. The default
                     setting is optimized for BWA. [50]
  -v, --version      show program's version number and exit
```

QTL-seq can be run from FASTQ (without or with trimming) and BAM. If you want to run QTL-seq from VCF, please use QTL-plot (example 5). Once you run QTL-seq, QTL-seq automatically completes the subprocesses.

+ [Example 1 : run QTL-seq from FASTQ without trimming](#example-1--run-qtl-seq-from-fastq-without-trimming)
+ [Example 2 : run QTL-seq from FASTQ with trimming](#example-2--run-qtl-seq-from-fastq-with-trimming)
+ [Example 3 : run QTL-seq from BAM](#example-3--run-qtl-seq-from-bam)
+ [Example 4 : run QTL-seq from multiple FASTQs and BAMs](#example-4--run-qtl-seq-from-multiple-fastqs-and-bams)
+ [Example 5 : run QTL-plot from VCF](#example-5--run-qtl-plot-from-vcf)



### Example 1 : run QTL-seq from FASTQ without trimming
```
qtlseq -r reference.fasta \
       -p parent.1.fastq.gz,parent.2.fastq.gz \
       -b1 bulk_1.1.fastq.gz,bulk_1.2.fastq.gz \
       -b2 bulk_2.1.fastq.gz,bulk_2.2.fastq.gz \
       -n1 20 \
       -n2 20 \
       -o example_dir
```

`-r` : reference fasta

`-p` : FASTQs of parent. Please input paired-end reads separated by a comma. FASTQ files can be gzipped.

`-b1` : FASTQs of bulk 1. Please input paired-end reads separated by a comma. FASTQ files can be gzipped.

`-b2` : FASTQs of bulk 2. Please input paired-end reads separated by a comma. FASTQ files can be gzipped.

`-n1` : number of individuals in bulk 1.

`-n2` : number of individuals in bulk 2.

`-o` : name of output directory. The specified directory should not already exist.

### Example 2 : run QTL-seq from FASTQ with trimming
```
qtlseq -r reference.fasta \
       -p parent.1.fastq.gz,parent.2.fastq.gz \
       -b1 bulk_1.1.fastq.gz,bulk_1.2.fastq.gz \
       -b2 bulk_2.1.fastq.gz,bulk_2.2.fastq.gz \
       -n1 20 \
       -n2 20 \
       -o example_dir \
       -T \
       -a adapter.fasta
```

`-r` : reference fasta

`-p` : FASTQs of parent. Please input paired-end reads separated by a comma. FASTQ files can be gzipped.

`-b1` : FASTQs of bulk 1. Please input paired-end reads separated by a comma. FASTQ files can be gzipped.

`-b2` : FASTQs of bulk 1. Please input paired-end reads separated by a comma. FASTQ files can be gzipped.

`-n1` : number of individuals in bulk 1.

`-n2` : number of individuals in bulk 2.

`-o` : name of output directory. The specified directory should not already exist.

`-T` : trim your reads using trimmomatic.

`-a` : FASTA of adapter sequences for trimmomatic.

If you are using TrueSeq3, you can find the adapter sequences in the [Github page of Trimmomatic](https://github.com/usadellab/Trimmomatic/blob/main/adapters/TruSeq3-PE-2.fa). [This thread](https://github.com/usadellab/Trimmomatic/issues/20) is also helpful to preprare the adapter file.

### Example 3 : run QTL-seq from BAM
```
qtlseq -r reference.fasta \
       -p parent.bam \
       -b1 bulk_1.bam \
       -b2 bulk_2.bam \
       -n1 20 \
       -n2 20 \
       -o example_dir
```

`-r` : reference fasta

`-p` : BAM of parent.

`-b1` : BAM of bulk 1.

`-b2` : BAM of bulk 2.

`-n1` : number of individuals in bulk 1.

`-n2` : number of individuals in bulk 2.

`-o` : name of output directory. The specified directory should not already exist.

### Example 4 : run QTL-seq from multiple FASTQs and BAMs
```
qtlseq -r reference.fasta \
       -p parent_1.1.fastq.gz,parent_1.2.fastq.gz \
       -p parent_1.bam \
       -b1 bulk_11.1.fastq.gz,bulk_11.2.fastq.gz \
       -b1 bulk_12.bam \
       -b1 bulk_13.bam \
       -b2 bulk_21.1.fastq.gz,bulk_21.2.fastq.gz \
       -b2 bulk_22.bam \
       -b2 bulk_23.bam \
       -n1 20 \
       -n2 20 \
       -o example_dir
```

QTL-seq automatically merges multiple FASTQ and BAM files. Of course, you can merge FASTQ or BAM files using `cat` or `samtools merge` before inputting them into QTL-seq. If you specify `-p` multiple times, please make sure that those files include only one individual. On the other hand, `-b1` and `-b2` can include more than one individuals because those are bulked samples. QTL-seq automatically classifies FASTQ and BAM files based on whether comma exits or not.

### Example 5 : run QTL-plot from VCF
```
qtlplot -h

usage: qtlplot -v <VCF> -n1 <INT> -n2 <INT> -o <OUT_DIR>
               [-F <INT>] [-t <INT>] [-w <INT>] [-s <INT>] [-D <INT>]
               [-d <INT>] [-N <INT>] [-m <FLOAT>] [-S <INT>] [-e <DATABASE>]
               [--igv] [--indel]

QTL-plot version 2.2.7

options:
  -h, --help            show this help message and exit
  -v , --vcf            VCF file which contains parent, bulk1, and bulk2
                        in this order. This VCF file must have AD field.
  -n1 , --N-bulk1       Number of individuals in bulk 1.
  -n2 , --N-bulk2       Number of individuals in bulk 2.
  -o , --out            Output directory. The specified directory can already
                        exist.
  -F , --filial         Filial generation. This parameter must be
                        more than 1. [2]
  -t , --threads        Number of threads. If a value less than 1 is specified,
                        QTL-seq will use the maximum available threads. [2]
  -w , --window         Window size (kb). [2000]
  -s , --step           Step size (kb). [100]
  -D , --max-depth      Maximum depth of variants to be used. [250]
  -d , --min-depth      Minimum depth of variants to be used. [8]
  -N , --N-rep          Number of replicates for simulations to generate
                        null distribution. [5000]
  -m , --min-SNPindex   Cutoff of minimum SNP-index for clear results. [0.3]
  -S , --strand-bias    Filter out spurious homozygous genotypes in the cultivar
                        based on strand bias. If ADF (or ADR) is higher than
                        this cutoff when ADR (or ADF) is 0, that SNP will be
                        filtered out. If you want to disable this filtering,
                        set this cutoff to 0. [7]
  -e , --snpEff         Predict causal variants using SnpEff. Check
                        available databases in SnpEff.
  --igv                 Output IGV format file to check results on IGV.
  --indel               Plot SNP-index with INDEL.
  --line-colors         Colors for threshold lines in plots. Specify a
                        comma-separated list in the order of SNP-index,
                        p95, and p99. ["#FE5F55,#6FD08C,#E3B505"]
  --dot-colors          Colors for dots in plots. Specify a
                        comma-separated list in the order of bulk1,
                        bulk2, and delta. ["#74D3AE,#FFBE0B,#40476D"]
  --fig-width           Width allocated in chromosome figure. [7.5]
  --fig-height          Height allocated in chromosome figure. [4.0]
  --white-space         White space between figures. (This option
                        only affects vertical direction.) [0.6]
  -f , --format         Specify the format of an output image.
                        eps/jpeg/jpg/pdf/pgf/png/rgba/svg/svgz/tif/tiff
  --version             show program's version number and exit
```
QTL-plot is included in QTL-seq. QTL-seq run QTL-plot after making VCF. Then, QTL-plot will work with default parameters. If you want to change some parameters, you can use VCF inside of `(OUT_DIR/30_vcf/QTL-seq.vcf.gz)` to retry plotting process like below.

```
qtlplot -v OUT_DIR/30_vcf/QTL-seq.vcf.gz \
        -o ANOTHER_DIR_NAME \
        -n1 20 \
        -n2 20 \
        -w 2000 \
        -s 100
```

#### Using QTL-plot for a VCF which was made by yourself
In this case:
1. Ensure that your VCF includes the AD format.
2. Ensure that your VCF includes columns for the parent, bulk1, and bulk2 in this order

If you got an error, please try to run QTL-seq from FASTQ or BAM before asking in issues.

## Outputs
Inside of `OUT_DIR` is like below.
```
├── 10_ref
│  ├── reference.fasta
│  ├── reference.fasta.amb
│  ├── reference.fasta.ann
│  ├── reference.fasta.bwt
│  ├── reference.fasta.fai
│  ├── reference.fasta.pac
│  └── reference.fasta.sa
├── 20_bam
│  ├── bulk1.filt.bam
│  ├── bulk1.filt.bam.bai
│  ├── bulk2.filt.bam
│  ├── bulk2.filt.bam.bai
│  ├── parent.filt.bam
│  └── parent.filt.bam.bai
├── 30_vcf
│  ├── qtlseq.vcf.gz
│  └── qtlseq.vcf.gz.tbi
├── 40_qtlseq
│  ├── bulk1_SNPindex.png
│  ├── bulk2_SNPindex.png
│  ├── delta_SNPindex.png
│  ├── sliding_window.tsv
│  ├── sliding_window.p95.tsv
│  ├── sliding_window.p99.tsv
│  ├── snp_index.tsv
│  ├── snp_index.p95.tsv
│  └── snp_index.p99.tsv
└── log
   ├── alignment.log
   ├── bcftools.log
   ├── bwa.log
   ├── samtools.log
   └── tabix.log
```
- If you run QTL-seq with trimming, you will get the directory of `00_fastq` which includes FASTQs after trimming.
- You can check the results in `40_QTL-seq`.
  + `snp_index.tsv` : columns in this order.
    - **CHROM** : chromosome name
    - **POSI** : position in chromosome
    - **VARIANT** : SNP or INDEL
    - **DEPTH 1** : depth of bulk 1
    - **DEPTH 2** : depth of bulk 2
    - **p99** : 99% confidence interval of simulated delta SNP-index (absolute value)
    - **p95** : 95% confidence interval of simulated delta SNP-index (absolute value)
    - **SNP-index 1** : real SNP-index of bulk 1
    - **SNP-index 2** : real SNP-index of bulk 2
    - **DELTA SNP-index** : real delta SNP-index (bulk2 - bulk1)
  + `sliding_window.tsv` : columns in this order.
    - **CHROM** : chromosome name
    - **POSI** : central position of window
    - **MEAN p99** : mean of p99 (absolute value)
    - **MEAN p95** : mean of p95 (absolute value)
    - **MEAN SNP-index 1** : mean SNP-index of bulk 1
    - **MEAN SNP-index 2** : mean SNP-index of bulk 2
    - **MEAN DELTA SNP-index** : mean delta SNP-index
  + `QTL-seq_plot.png` : resulting plot (like below)
    - **<span style="color: blue; ">BLUE dot</span>** : variant
    - **<span style="color: red; ">RED line</span>** : mean SNP-index
    - **<span style="color: orange; ">ORANGE line</span>** : mean p99
    - **<span style="color: green; ">GREEN line</span>** : mean p95

<img src="https://github.com/YuSugihara/QTL-seq/blob/master/images/2_result.png" width=600>

## Additional Resources

### QTL-seq commands breakdown

For a detailed breakdown of the commands used in QTL-seq, including explanations of each step, parameters, and best practices for troubleshooting, please refer to the [QTL-seq Commands Breakdown](docs/01_qtlseq_commands_breakdown.md) document.

## About multiple testing correction
This function has been deprecated since v2.2.5.
We highly recommend running QTL-seq without this function.
However, if you would like to use this function, you can use it with versions of QTL-seq older than v2.2.5.

## Build and use your own database for snpEff
If you want to use your own database for snpEff, you need additional steps.
Here we assume that you installed QTL-seq via anaconda distribution, creating new environment with `conda create`.

1. Find the directory of snpEff that includes snpEff script, configuration file and database. You can find it in `/home/anaconda3/envs/{your_env_name_installed_qtlseq}/share/snpeff-5.0-0/`. `anaconda3` may be `miniconda3`. Also, the version of snpeff may be different.

2. Go to this directory and follow the snpEff manual to build the database.
Don't forget to add your database info to the snpEff configuration file.
https://pcingola.github.io/SnpEff/se_buildingdb/#add-a-genome-to-the-configuration-file

3. Run QTL-seq with option `-e {your_database_name}`
