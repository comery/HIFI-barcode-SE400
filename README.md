# HIFI-barcode-SE400
The BGISEQ-500 platform has launched a new test sequencing kits capable of single-end 400 bp sequencing (SE400), which offers a simple and reliable way to achieve DNA barcodes efficiently. In this study, we explored the potential of the BGISEQ-500 SE400 sequencing in DNA barcode reference construction, meanwhile provided an updated HIFI-Barcode software package that can generate COI barcode assemblies using HTS reads of length >= 400 bp.


### Manual
Download manual book to see details. => [manual book](https://github.com/comery/HIFI-barcode-SE400/raw/master/HIFI-SE_manual.pdf)

### Changelog

#### latest version 2.0.2 Python3
- v2.0.2 2020-06-25 add a function in "taxonomy" module that can BLAST the 5’ and 3’ end of the barcode sequences and then compare taxonomies for consistency to further validate the assembly accuracy.
- v2.0.1 2019-04-15 fix a bug of "assign", and degenerated base in primer now is accepted.
- v1.0.5 2019-04-09 add support to compressed fastq, fix a bug of taxonomy
- v1.0.4 2019-04-02 Fix a bug of "polish", and update the bold_identification module
- v1.0.3 2018-12-14 Fix a bug of "trim"
- v1.0.2 2018-12-10 Add "-trim" function in filter;
	accept mismatches in tag or primer sequence,
   	when demultiplexing; accept uneven reads to
   	assembly; add "-ds" to drop short reads before
   	assembly.
- v1.0.1 2018-12-2  Add "polish" function
- v1.0.0
	HIFI-SE v1.0.0 2018/11/22. Changers form previous version:
	- Formatted python code writing style as PEP8.
	- Fixed several small bugs.
- v0.0.3
	HIFI-SE v0.03 2018/11/15. Changes from previous version:
	- Modify the description of some arguments for better understanding.
- v0.0.1  
	HIFI-SE v0.0.1 2018/11/03 beat version, establish the framework and archive almost complete functions.

#### original Perl version & Python, rough sources
>0.expected_error.pl   
1.split_extract.pl  
2.hificonnect.pl

>0.expected_error.py  
1.split_extract.py  
2.hificonnect.py  

### Installation

#### System requirement and dependencies 
Operating system: HIFI-SE is designed to run on most platforms, including UNIX, Linux and MacOS/X. Microsoft Windows. We have tested on Linux and on MacOS/X, because these are the machines we develop on. HIFI-SE is written in python language, and a version 3.5 or higher is required.
#### Dependencies: 
- biopython version 1.5 or higher (required). Please check https://biopython.org/ and https://pypi.org/project/biopython/#description for more details on installation of biopython.
- Another python package - bold_identification is also required for getting complete function of HIFI-SE. See https://pypi.org/project/bold-identification/
- HIFI-SE supposed you have installed the VSEARCH on your device, and its path in your $PATH. See https://github.com/torognes/vsearch

#### Install

1. [Recommended ] I only deploy my latest version on github, so you can clone this repository to your local computer. However, it would not solve package dependencies, thus you need to install biopython and bold_identification before using HIFI-SE software.<b>(NOTE: pip is a link from pip3)</b>

	```shell
	git clone https://github.com/comery/HIFI-barcode-SE400.git
	pip install biopython
	pip install bold_identification  
	```
	
2. [Perhaps out of date] Installation by pip is recommended because it will solve package dependencies automatically, including biopython and bold_identification packages. 

	```pip install HIFI-SE``` 


### Usage (latest)

```shell
python3 HIFI-SE.py
```
or 

```shell
HIFI-SE.py
```

```text
usage: HIFI-SE [-h] [-v]
               {all,filter,assign,assembly,polish,bold_identification} ...

Description

    An automatic pipeline for HIFI-SE400 project, including filtering
    raw reads, assigning reads to samples, assembly HIFI barcodes
    (COI sequences), polished assemblies, and do tax identification.
    See more: https://github.com/comery/HIFI-barcode-SE400

Versions

    2.0.2 (20200625)

Authors

    yangchentao at genomics.cn, BGI.
    mengguanliang at genomics.cn, BGI.

positional arguments:
  {all,filter,assign,assembly,polish,bold_identification}
    all                 run filter, assign and assembly.
    filter              remove or trim reads with low quality.
    assign              assign reads to samples by tags.
    assembly            do assembly from assigned reads,
                        output raw HIFI barcodes.
    polish              polish COI barcode assemblies,
                        output confident barcodes.
    taxonomy            do taxa identification on BOLD system

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

#### run by steps [filter -> assign -> assembly]

- ```python3 HIFI-SE.py filter ```

```text
usage: HIFI-SE filter [-h] -outpre <STR> -raw <STR> [-phred <INT>] [-e <INT>]
                      [-q <INT> <INT>] [-trim] [-n <INT>]

optional arguments:
  -h, --help      show this help message and exit

common arguments:
  -outpre <STR>   prefix for output files

filter arguments:
  -raw <STR>      input raw Single-End fastq file, and only
                  adapters should be removed; supposed on
                  Phred33 score system (BGISEQ-500)
  -phred <INT>    Phred score system, 33 or 64, default=33
  -e <INT>        expected error threshod, default=10
                  see more: http://drive5.com/usearch/manual/exp_errs.html
  -q <INT> <INT>  filter by base quality; for example: '20 5' means
                  dropping read which contains more than 5 percent of
                  quality score < 20 bases.
  -trim           whether to trim 5' end of read, it adapts to -e mode
                  or -q mode
  -n <INT>        remove reads containing [INT] Ns, default=1
```

- ```python3 HIFI-SE.py assign```

```text
usage: HIFI-SE assign [-h] -outpre <STR> -index INT -fq <STR> -primer <STR>
                      [-outdir <STR>] [-tmis <INT>] [-pmis <INT>]

optional arguments:
  -h, --help     show this help message and exit

common arguments:
  -outpre <STR>  prefix for output files

index arguments:
  -index INT     the length of tag sequence in the ends of primers

when only run assign arguments:
  -fq <STR>      cleaned fastq file (*.fq.gz, *.fq)

assign arguments:
  -primer <STR>  taged-primer list, on following format:
                 Rev001   AAGCTAAACTTCAGGGTGACCAAAAAATCA
                 For001   AAGCGGTCAACAAATCATAAAGATATTGG
                 ...
                 this format is necessary, degenerated base is also accepted.
  -outdir <STR>  output directory for assignment,default="assigned"
  -tmis <INT>    mismatch tolerance in tag when demultiplexing, default=0
  -pmis <INT>    mismatch tolerance in primer when demultiplexing, default=1
```
- ```python3 HIFI-SE.py assembly```

```
usage: HIFI-SE assembly [-h] -outpre <STR> -index INT -list FILE
                        [-vsearch <STR>] [-threads <INT>] [-cid FLOAT]
                        [-min INT] [-max INT] [-oid FLOAT] [-tp INT] [-ab INT]
                        [-seqs_lim INT] [-len INT] [-ds] [-mode INT] [-rc]
                        [-codon INT] [-frame INT]

optional arguments:
  -h, --help      show this help message and exit

common arguments:
  -outpre <STR>   prefix for output files

index arguments:
  -index INT      the length of tag sequence in the ends of primers

only run assembly arguments(not all):
  -list FILE      input file, fastq file list. [required]

software path:
  -vsearch <STR>  vsearch path(only needed if vsearch is not in $PATH)
  -threads <INT>  threads for vsearch, default=2
  -cid FLOAT      identity for clustering, default=0.98

assembly arguments:
  -min INT        minimun length of overlap, default=80
  -max INT        maximum length of overlap, default=90
  -oid FLOAT      minimun similarity of overlap region, default=0.95
  -tp INT         how many clusters will be used inassembly, recommend 2
  -ab INT         keep clusters to assembly if its abundance >=INT
  -seqs_lim INT   reads number limitation. by default,
                  no limitation for input reads
  -len INT        standard read length, default=400
  -ds             drop short reads away before assembly
  -mode INT       1 or 2; modle 1 is to cluster and keep
                  most [-tp] abundance clusters, or clusters
                  abundance more than [-ab], and then make a
                  consensus sequence for each cluster.
                  modle 2 is directly to make only one consensus
                  sequence without clustering. default=1
  -rc             whether to check amino acid
                  translation for reads, default not

translation arguments(when set -rc or -cc):
  -codon INT      codon usage table used to checktranslation, default=5
  -frame INT      start codon shift for amino acidtranslation, default=1
```

- ```python3 HIFI-SE.py taxonomy```

```text
usage: taxonomy [-h] -i <str> [-f <str>] -o <str>
                [-d {COX1,COX1_SPECIES,COX1_SPECIES_PUBLIC,COX1_L640bp,ITS,MATK_RBCL}]
                [-n <int>] [-r <int>] [-c] [-D] [--version]

To identify taxa of given sequences from BOLD (http://www.boldsystems.org/).
Some sequences can fail to get taxon information, which can be caused by
TimeoutException if your network to the BOLD server is bad.
Those sequences will be output in the file '*.TimeoutException.fasta'.

You can:
1) run another searching with the same command directly (but add -c option);
2) lengthen the time to wait for each query (-t option);
3) increase submission times (-r option) for a sequence.

Also, the sequences without BOLD matches will be output in the
file '*.NoBoldMatchError.fasta'

By mengguanliang AT genomics DOT cn.
See https://github.com/linzhi2013/bold_identification.

version: 0.0.25

optional arguments:
  -h, --help            show this help message and exit
  -i <str>              input file name
  -f <str>              input file format [fasta]
  -o <str>              outfile
  -d {COX1,COX1_SPECIES,COX1_SPECIES_PUBLIC,COX1_L640bp,ITS,MATK_RBCL}
                        database to search [COX1]
  -n <int>              how many first top hits will be output. [1]
  -r <int>              Maximum submission time for a sequence, useful for
                        handling TimeOutException. [4]
  -c                    continuous mode, jump over the ones already in "-o"
                        file, will resubmit all the remained. use "-cc" to
                        also jump over the ones in "*.NoBoldMatchError.fasta"
                        file. [0]
  -D                    debug mode output [False]
  --version             show program's version number and exit
  
```


### Quickstart
#### Files used in tutorial
All related files could be found from here. The important files for tutorial are:

- <i>raw.fastq.gz</i>, raw output fastq file generated from BGISEQ-500 SE400 module.
- <i>indexed_primer.list</i>, tagged primer list


#### Run in "all"
Example:

```shell
python3 HIFI-SE.py all -outpre hifi -trim -e 5 -raw test.raw.fastq -index 5 -primer index_primer.list -mode 1 -cid 0.98 -oid 0.95 -seqs_lim 50000 -threads 4 -tp 2 
```

### Citation

Current:

Yang, C., Zheng, Y., Tan, S., Meng, G., Rao, W., Yang, C., ... & Liu, S. (2020). Efficient COI barcoding using high throughput single-end 400 bp sequencing. BMC genomics, 21(1), 1-10.

PrePrint:

Yang, Chentao, Shangjin Tan, Guanliang Meng, David G. Bourne, Paul A. O'brien, Junqiang Xu, Sha Liao, Ao Chen, Xiaowei Chen, and Shanlin Liu. "Access COI barcode efficiently using high throughput Single End 400 bp sequencing." BioRxiv (2018): 498618.



