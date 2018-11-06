# HIFI-barcode-SE400
The BGISEQ-500 platform has launched a new test sequencing kits capable of single-end 400 bp sequencing (SE400), which offers a simple and reliable way to achieve DNA barcodes efficiently. In this study, we explored the potential of the BGISEQ-500 SE400 sequencing in DNA barcode reference construction, meanwhile provided an updated HIFI-Barcode software package that can generate COI barcode assemblies using HTS reads of length > 400 bp.


### Versions
#### version 1.2 Perl & Python
>0.expected_error.pl   
1.split_extract.pl  
2.hificonnect.pl

>0.expected_error.py  
1.split_extract.py  
2.hificonnect.py  

#### version 2.0 Python
>hifi.v2.py

### Usage (latest)

```shell
python3 version2.0/hifi.v2.py
```

```text
usage: HIFI-SE400.py [-h] {all,filter,assign,assembly} ...

Description An automatic pipeline for HIFI-SE400 project, including filtering
raw reads, assigning reads to samples, assembly HIFI barcodes (COI sequences).
Version 2.0 Author yangchentao at genomics.cn, BGI.

positional arguments:
  {all,filter,assign,assembly}
    all                 run filter, assign and assembly
    filter              filter raw reads
    assign              assign reads to samples
    assembly            do assembly from input fastq reads, output HIFI
                        barcodes.

optional arguments:
  -h, --help            show this help message and exit

```

#### run in "all"
Example:

```shell
python3 version2.0/hifi.v2.py all -outpre hifi -raw test.raw.fastq -index 5 -primer index_primer.list -cid 0.98 -oid 0.95 -seqs_lim 50000 -threads 4 -tp 2
```
#### run by steps [filter -> assign -> assembly]

- ```python3 version2.0/hifi.v2.py filter ```

```text
usage: HIFI-SE400.py filter [-h] -outpre <STR> -raw <STR> [-e <INT>]

optional arguments:
  -h, --help     show this help message and exit

common arguments:
  -outpre <STR>  outprefix for process.

Filter arguments:
  -raw <STR>     raw fastq file
  -e <INT>       expected error number threshod [10]
```

- ```python3 hifi.v2.py assign```

```text
uusage: HIFI-SE400.py assign [-h] -outpre <STR> -index INT -fq <STR> -primer
                            <STR> [-outdir <STR>]

optional arguments:
  -h, --help     show this help message and exit

common arguments:
  -outpre <STR>  outprefix for process.

index arguments:
  -index INT     index sequence lenght

when only run assign arguments:
  -fq <STR>      cleaned fastq file

assign arguments:
  -primer <STR>  taged primer list, like following lines: Rev001
                 AAGCTAAACTTCAGGGTGACCAAAAAATCAFor001
                 AAGCGGTCAACAAATCATAAAGATATTGG...this format is necessary!
  -outdir <STR>  output directory for assignment
```
- ```python3 hifi.v2.py assembly```

```
usage: HIFI-SE400.py assembly [-h] -outpre <STR> -index INT -list FILE
                              [-vsearch <STR>] [-threads <INT>] [-cid FLOAT]
                              [-min INT] [-max INT] [-oid FLOAT] [-tp INT]
                              [-ab INT] [-seqs_lim INT] [-len INT] [-mode INT]
                              [-rc] [-cc] [-codon INT] [-frame INT]

optional arguments:
  -h, --help      show this help message and exit

common arguments:
  -outpre <STR>   outprefix for process.

index arguments:
  -index INT      index sequence lenght

when only run assembly arguments:
  -list FILE      input file, fastq file list. [required]

software path:
  -vsearch <STR>  vsearch path directory (only needed if vsearch is not in
                  PATH)
  -threads <INT>  threads for vsearch
  -cid FLOAT      clustering identity rate [0.98]

assembly arguments:
  -min INT        minimun length of overlap [80]
  -max INT        maximum length of overlap [90]
  -oid FLOAT      cutoff of identity of overlap region [0.95]
  -tp INT         how many clusters using in assembly.
  -ab INT         keep all clusters to assembly if abundance >=INT
  -seqs_lim INT   reads number limitation. [0]
  -len INT        standard reads length [400]
  -mode INT       modle 1 is to cluster and keep 3 clusters with most
                  abundancemodle 2 is to make a consensus sequence using all
                  reads or allcodon checked reads if you set "-rc" * --cid is
                  invaild in this mode
  -rc             whether to check reads' codon translation
  -cc             whether to check final COI contig's condon translation
  -codon INT      codon table using to check translation [5], by the way,
                  table [4,5] have same effect for COI gene.
  -frame INT      translation start shift [1]
```


