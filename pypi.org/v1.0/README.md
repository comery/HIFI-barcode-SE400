# HIFI-barcode-SE400
The BGISEQ-500 platform has launched a new test sequencing kits capable of single-end 400 bp sequencing (SE400), which offers a simple and reliable way to achieve DNA barcodes efficiently. In this study, we explored the potential of the BGISEQ-500 SE400 sequencing in DNA barcode reference construction, meanwhile provided an updated HIFI-Barcode software package that can generate COI barcode assemblies using HTS reads of length > 400 bp.


### Versions

new release: 1.0 2018/11/13

### Usage (latest)


```shell
HIFI-SE
```

```text
usage: HIFI-SE [-h] {all,filter,assign,assembly,bold_identification} ...

Description

	An automatic pipeline for HIFI-SE400 project, including filtering raw reads,
	assigning reads to samples, assembly HIFI barcodes (COI sequences).

Version
	1.0 2018-11-3

Author
	yangchentao at genomics.cn, BGI.
	mengguanliang at genomics.cn, BGI.


positional arguments:
  {all,filter,assign,assembly,bold_identification}
    all                 run filter, assign and assembly
    filter              filter raw reads
    assign              assign reads to samples
    assembly            do assembly from input fastq reads,
                        output HIFI barcodes.
    bold_identification
                        do taxa identification on BOLD system,

optional arguments:
  -h, --help            show this help message and exit

```

#### run in "all"
Example:

```shell
HIFI-SE all -outpre hifi -raw test.raw.fastq -index 5 -primer index_primer.list -cid 0.98 -oid 0.95 -seqs_lim 50000 -threads 4 -tp 2
```
#### run by steps [filter -> assign -> assembly]

- ```python3 HIFI-SE.py filter ```

```text
usage: HIFI-SE filter [-h] -outpre <STR> -raw <STR> [-e <INT>]
                      [-q <INT> <INT>] [-n <INT>]

optional arguments:
  -h, --help      show this help message and exit

common arguments:
  -outpre <STR>   outprefix for process.

filter arguments:
  -raw <STR>      input raw singled-end fastq file, (Phred33)
  -e <INT>        expected error number threshod, P = 10–Q/10, default=10
  -q <INT> <INT>  filter by quality method,  Q = –10 log10(P),
                  filter out low quality reads. example: 20 5, it means
                  dropping read which contains more than 5 percent of
                  quality score < 20 bases.
  -n <INT>        remove reads containing [INT] Ns, default=1
```

- ```python3 HIFI-SE.py assign```

```text
usage: HIFI-SE assign [-h] -outpre <STR> -index INT -fq <STR> -primer <STR>
                      [-outdir <STR>]

optional arguments:
  -h, --help     show this help message and exit

common arguments:
  -outpre <STR>  outprefix for process.

index arguments:
  -index INT     index sequence lenght

when only run assign arguments:
  -fq <STR>      cleaned fastq file

assign arguments:
  -primer <STR>  taged primer list, like following lines:
                 Rev001	AAGCTAAACTTCAGGGTGACCAAAAAATCA
                 For001	AAGCGGTCAACAAATCATAAAGATATTGG
                 ...
                 this format is necessary!
  -outdir <STR>  output directory for assignment
```
- ```python3 HIFI-SE.py assembly```

```
usage: HIFI-SE assembly [-h] -outpre <STR> -index INT -list FILE
                        [-vsearch <STR>] [-threads <INT>] [-cid FLOAT]
                        [-min INT] [-max INT] [-oid FLOAT] [-tp INT] [-ab INT]
                        [-seqs_lim INT] [-len INT] [-mode INT] [-rc] [-cc]
                        [-codon INT] [-frame INT]

optional arguments:
  -h, --help      show this help message and exit

common arguments:
  -outpre <STR>   outprefix for process.

index arguments:
  -index INT      index sequence lenght

when only run assembly arguments:
  -list FILE      input file, fastq file list. [required]

software path:
  -vsearch <STR>  vsearch path(only needed if vsearch is not in $PATH)
  -threads <INT>  threads for vsearch
  -cid FLOAT      identity for clustering [0.98]

assembly arguments:
  -min INT        minimun length of overlap [80]
  -max INT        maximum length of overlap [90]
  -oid FLOAT      minimun identity of overlap region [0.95]
  -tp INT         how many clusters using in assembly. default=2
  -ab INT         keep all clusters to assembly if its abundance >=INT
  -seqs_lim INT   reads number limitation. [0]
  -len INT        standard reads length [400]
  -mode INT       modle 1 is to cluster and keep most [-tp] abundance clusters,
                  or clusters abundance more than [-ab], and then make a consensus
                  sequence for each cluster. modle 2 is directly to make only
                  consensus sequence without clustering.
  -rc             whether to check amino acid translation for reads
  -cc             whether to check final COI contig's amino acid translation
  -codon INT      codon table using to check translation [5],
                  by the way, table [4,5] have same effect for COI gene.
  -frame INT      translation start shift [1]
```

#### Github page
https://github.com/comery/HIFI-barcode-SE400
