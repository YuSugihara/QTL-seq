# QTL-seq User Guide
#### version 2.0.8

## Table of contents
- [What is QTL-seq?](#What-is-QTL-seq)
- [Installation](#Installation)
  + [Dependencies](#Dependencies)
  + [Installation using bioconda](#Installation-using-bioconda)
  + [Manual Installation](#Manual-Installation)
- [Usage](#Usage)
  + [Example 1 : run QTL-seq from FASTQ without trimming](#Example-1--run-QTL-seq-from-FASTQ-without-trimming)
  + [Example 2 : run QTL-seq from FASTQ with trimming](#Example-2--run-QTL-seq-from-FASTQ-with-trimming)
  + [Example 3 : run QTL-seq from BAM](#Example-3--run-QTL-seq-from-BAM)
  + [Example 4 : run QTL-seq from multiple FASTQs and BAMs](#Example-4--run-QTL-seq-from-multiple-FASTQs-and-BAMs)
  + [Example 5 : run QTL-plot from VCF](#Example-5--run-QTL-plot-from-VCF)
- [Outputs](#Outputs)

## What is QTL-seq?
<img src="https://210a94ef-a-fc7f2be1-s-sites.googlegroups.com/a/ibrc.or.jp/genome-e/home/bioinformatics-team/mutmap/QTL-seq_LOGO.jpg?attachauth=ANoY7coJoWWEVItHoU2_q2dDtQkdBtsHOiv228_rPBP8ZvISV9vvS7z-YRmxG5_JLCMhwPrMg92SBMU-ULQd5F1X2Tioz5Vvv_gLecE3WSMlrBYu5y69PHEuhTWWDM8hRruMJWCrVIdjeQoDL26KWbZUFVEZOHHXX5zL-s2B0UqH3zKDxurCfpmrg_gbE7y_8D9gvaGEAYe73HOR1Jl7WjpdeYijWeqanQzUnLwWzMnITpKxLzXD7fD5ebjBLI5wNZo2j7UTfRIv&attredirects=0" width=200>

Bulked segregant analysis, as implemented in  QTL-seq (Takagi et al., 2013), is a powerful and efficient method to identify agronomically important loci in crop plants. QTL-seq was adapted from MutMap to identify quantitative trait loci. It utilizes sequences pooled from two segregating progeny populations with extreme opposite traits (e.g. resistant vs susceptible) and a single whole-genome resequencing of either of the parental cultivars. While the original QTL-seq algorithm did not assume a highly heterozygous genome, a “modified QTL-seq” has been developed to handle this situation using high resolution mapping ([Itoh et al., 2019](https://doi.org/10.1007/s00122-019-03396-z)).

#### Citation
- Hiroki Takagi, Akira Abe, Kentaro Yoshida, Shunichi Kosugi, Satoshi Natsume, Chikako Mitsuoka, Aiko Uemura, Hiroe Utsushi, Muluneh Tamiru, Shohei Takuno, Hideki Innan, Liliana M. Cano, Sophien Kamoun, Ryohei Terauchi (2013).  QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations. Plant journal 74:174-183. [[URL]](https://doi.org/10.1111/tpj.12105)
- Yu Sugihara, Lester Young, Hiroki Yaegashi, Satoshi Natsume, Daniel J. Shea, Hiroki Takagi, Helen Booker, Ryohei Terauchi, Akira Abe (in preparation). High performance pipeline for MutMap and QTL-seq.

## Installation
### Dependencies
#### Softwares
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

### Installation using bioconda
You can install QTL-seq using [bioconda](https://bioconda.github.io/index.html).
```
$ conda install -c bioconda qtlseq
```

### Mannual Installation
If you got a error during installation, you can install QTL-seq, manually.
```
$ git clone https://github.com/YuSugihara/QTL-seq.git
$ cd QTL-seq
$ pip install -e .
```
Then you have to install other dependencies by yourself. We highly recommend you to install SnpEff and Trimmomatic using bioconda.
```
$ conda install -c bioconda snpeff
$ conda install -c bioconda trimmomatic
```
After installation, please check whether SnpEff and Trimmomatic work through the commands like below.
```
$ snpEff --help
$ trimmomatic --help
```

## Usage
```
$ qtlseq -h

usage: qtlseq -r <FASTA> -p <BAM|FASTQ> -b1 <BAM|FASTQ>
              -b2 <BAM|FASTQ> -n1 <INT> -n2 <INT> -o <OUT_DIR>
              [-F <INT>] [-T] [-e <DATABASE>]

QTL-seq version 2.0.8

optional arguments:
  -h, --help         show this help message and exit
  -r , --ref         Reference fasta.
  -p , --parent      fastq or bam of parent. If you specify
                     fastq, please separate pairs by comma,
                     e.g. -p fastq1,fastq2. You can use this
                     optiion multiple times
  -b1 , --bulk1      fastq or bam of bulk 1. If you specify
                     fastq, please separate pairs by comma,
                     e.g. -b1 fastq1,fastq2. You can use this
                     optiion multiple times
  -b2 , --bulk2      fastq or bam of bulk 2. If you specify
                     fastq, please separate pairs by comma,
                     e.g. -b2 fastq1,fastq2. You can use this
                     optiion multiple times
  -n1 , --N-bulk1    Number of individuals in bulk 1.
  -n2 , --N-bulk2    Number of individuals in bulk 2.
  -o , --out         Output directory. Specified name must not
                     exist.
  -F , --filial      Filial generation. This parameter must be
                     more than 1. [2]
  -t , --threads     Number of threads. If you specify the number
                     below one, then QTL-seq will use the threads
                     as many as possible. [2]
  -w , --window      Window size (kb). [2000]
  -s , --step        Step size (kb). [100]
  -D , --max-depth   Maximum depth of variants which will be used. [250]
  -d , --min-depth   Minimum depth of variants which will be used. [8]
  -N , --N-rep       Number of replicates for simulation to make
                     null distribution. [5000]
  -T, --trim         Trim fastq using trimmomatic.
  -a , --adapter     FASTA of adapter sequences. This will be used
                     when you specify "-T" for trimming.
  --trim-params      Parameters for trimmomatic. Input parameters
                     must be separated by comma with following
                     order: phred, ILLUMINACLIP, LEADING, TRAILING,
                     SLIDINGWINDOW, MINLEN. If you want to remove
                     adapters of illumina, please specify FASTA of
                     the adapter sequences with "--adapter". Specified
                     name will be inserted into <ADAPTER_FASTA>. If you
                     don't specify it, adapter trimming will be skipped.
                     [33,<ADAPTER_FASTA>:2:30:10,20,20,4:15,75]
  -e , --snpEff      Predict causal variant using SnpEff. Please
                     check available databases in SnpEff.
  --mem              Maximum memory per thread when bam sorted;
                     suffix K/M/G recognized. [1G]
  -q , --min-MQ      Minimum mapping quality in mpileup. [40]
  -Q , --min-BQ      Minimum base quality in mpileup. [18]
  -C , --adjust-MQ   "adjust-MQ" in mpileup. Default parameter
                     is suited for BWA. [50]
  -v, --version      show program's version number and exit
```

QTL-seq can run from FASTQ (without or with trimming) and BAM. If you want to run QTL-seq from VCF, please use QTL-plot (example 5). Once you run QTL-seq, QTL-seq automatically complete the subprocesses.

+ [Example 1 : run QTL-seq from FASTQ without trimming](#Example-1--run-QTL-seq-from-FASTQ-without-trimming)
+ [Example 2 : run QTL-seq from FASTQ with trimming](#Example-2--run-QTL-seq-from-FASTQ-with-trimming)
+ [Example 3 : run QTL-seq from BAM](#Example-3--run-QTL-seq-from-BAM)
+ [Example 4 : run QTL-seq from multiple FASTQs and BAMs](#Example-4--run-QTL-seq-from-multiple-FASTQs-and-BAMs)
+ [Example 5 : run QTL-plot from VCF](#Example-5--run-QTL-plot-from-VCF)



### Example 1 : run QTL-seq from FASTQ without trimming
```
$ qtlseq -r reference.fasta \
         -p parent.1.fastq,parent.2.fastq \
         -b1 bulk_1.1.fastq,bulk_1.2.fastq \
         -b2 bulk_2.1.fastq,bulk_2.2.fastq \
         -n1 20 \
         -n2 20 \
         -o example_dir
```

`-r` : reference fasta

`-p` : FASTQs of parent. Please input pair-end reads separated by comma. FASTQs can be gzipped.

`-b1` : FASTQs of bulk 1. Please input pair-end reads separated by comma. FASTQs can be gzipped.

`-b2` : FASTQs of bulk 2. Please input pair-end reads separated by comma. FASTQs can be gzipped.

`-n1` : number of individuals in bulk 1.

`-n2` : number of individuals in bulk 2.

`-o` : name of output directory. Specified name cannot exist.

### Example 2 : run QTL-seq from FASTQ with trimming
```
$ qtlseq -r reference.fasta \
         -p parent.1.fastq,parent.2.fastq \
         -b1 bulk_1.1.fastq,bulk_1.2.fastq \
         -b2 bulk_2.1.fastq,bulk_2.2.fastq \
         -n1 20 \
         -n2 20 \
         -o example_dir \
         -T
```

`-r` : reference fasta

`-p` : FASTQs of parent. Please input pair-end reads separated by comma. FASTQs can be gzipped.

`-b1` : FASTQs of bulk 1. Please input pair-end reads separated by comma. FASTQs can be gzipped.

`-b2` : FASTQs of bulk 1. Please input pair-end reads separated by comma. FASTQs can be gzipped.

`-n1` : number of individuals in bulk 1.

`-n2` : number of individuals in bulk 2.

`-o` : name of output directory. Specified name cannot exist.

`-T` : trim your reads using trimmomatic.

### Example 3 : run QTL-seq from BAM
```
$ qtlseq -r reference.fasta \
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

`-o` : name of output directory. Specified name cannot exist.

### Example 4 : run QTL-seq from multiple FASTQs and BAMs
```
$ qtlseq -r reference.fasta \
         -p parent_1.1.fastq,parent_1.2.fastq \
         -p parent_1.bam \
         -b1 bulk_11.1.fastq,bulk_11.2.fastq \
         -b1 bulk_12.bam \
         -b1 bulk_13.bam \
         -b2 bulk_21.1.fastq,bulk_21.2.fastq \
         -b2 bulk_22.bam \
         -b2 bulk_23.bam \
         -n1 20 \
         -n2 20 \
         -o example_dir
```

QTL-seq can automatically merge multiple FASTQs and BAMs. Of course, you can merge FASTQs or BAMs using `cat` or `samtools merge` before input them to QTL-seq. If you specify `-p` multiple times, please make sure that those files include only 1 individual. On the other hand, `-b1` and `-b2` can include more than 1 individuals because those are bulked samples. QTL-seq can automatically classify FASTQs and BAMs from whether comma exits or not.

### Example 5 : run QTL-plot from VCF
```
$ qtlplot -h

usage: qtlplot -v <VCF> -n1 <INT> -n2 <INT> -o <OUT_DIR>
               [-F <INT>] [-t <INT>] [-w <INT>] [-s <INT>] [-D <INT>]
               [-d <INT>] [-N <INT>] [-m <FLOAT>] [-S <INT>] [-e <DATABASE>]
               [--igv] [--corr] [--indel]

QTL-plot version 2.0.8

optional arguments:
  -h, --help            show this help message and exit
  -v , --vcf            VCF which contains parent, bulk1 and bulk2
                        in this order.
  -n1 , --N-bulk1       Number of individuals in bulk 1.
  -n2 , --N-bulk2       Number of individuals in bulk 2.
  -o , --out            Output directory. Specified name can exist.
  -F , --filial         Filial generation. This parameter must be
                        more than 1. [2]
  -t , --threads        Number of threads. If you specify the number
                        below one, then QTL-plot will use the threads
                        as many as possible. [2]
  -w , --window         Window size (kb). [2000]
  -s , --step           Step size (kb). [100]
  -D , --max-depth      Maximum depth of variants which will be used. [250]
  -d , --min-depth      Minimum depth of variants which will be used. [8]
  -N , --N-rep          Number of replicates for simulation to make
                        null distribution. [5000]
  -m , --min-SNPindex   Cutoff of minimum SNP-index for clear results. [0.3]
  -S , --strand-bias    Filter spurious homo genotypes in cultivar using
                        strand bias. If ADF (or ADR) is higher than this
                        cutoff when ADR (or ADF) is 0, that SNP will be
                        filtered out. If you want to supress this filtering,
                        please set this cutoff to 0. [7]
  -e , --snpEff         Predict causal variant using SnpEff. Please
                        check available databases in SnpEff.
  --igv                 Output IGV format file to check results on IGV.
  --corr                Use corrected threshold in Huang et al. (2019).
                        If you specify this option, the options related to
                        the simulation ("-N" and "-F") will be ignored.
  --indel               Plot SNP-index with INDEL.
  --fig-width           Width allocated in chromosome figure. [7.5]
  --fig-height          Height allocated in chromosome figure. [4.0]
  --white-space         White space between figures. (This option
                        only affect vertical direction.) [0.6]
  --version             show program's version number and exit
```
QTL-plot is included in QTL-seq. QTL-seq run QTL-plot after making VCF. Then, QTL-plot will work with default parameters. If you want to change some parameters, you can use VCF inside of `(OUT_DIR/30_vcf/QTL-seq.vcf.gz)` to retry plotting process like below.

```
$ qtlplot -v OUT_DIR/30_vcf/QTL-seq.vcf.gz \
          -o ANOTHER_DIR_NAME \
          -n1 20 \
          -n2 20 \
          -w 2000 \
          -s 100
```

#### Use QTL-plot for VCF which was made by yourself
In this case, please make sure that:
1. Your VCF include AD format.
2. Your VCF include three columns of parent, bulk1 and bulk2 in this order.

If you got a error, please try to run QTL-seq from FASTQ or BAM before asking in issues.

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
│  └── snp_index.tsv
└── log
   ├── bcftools.log
   ├── bgzip.log
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

<img src="https://user-images.githubusercontent.com/34593586/72580948-00e60500-3921-11ea-850f-dcf9a8d75e74.png" width=600>
