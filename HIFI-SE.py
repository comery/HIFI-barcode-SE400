#!/usr/bin/env python3
import os
import sys
import time
import gzip
import argparse
import subprocess
from bold_identification.BOLD_identification import (
    main as bold_identification,
)

t = time.time()

try:
    import Bio
except:
    sys.exit("package biopython not found! Please install it!")
else:
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna


###############################################################################
#####------------------------- parameters --------------------------------#####

## common group  ##
common_parser = argparse.ArgumentParser(add_help=False)

common_group = common_parser.add_argument_group("common arguments")

common_group.add_argument(
    "-outpre",
    metavar="<STR>",
    required=True,
    help="prefix for output files",
)
## index group ##
index_parser = argparse.ArgumentParser(add_help=False)

index_group = index_parser.add_argument_group("index arguments")

index_group.add_argument(
    "-index",
    metavar="INT",
    type=int,
    required=True,
    help="the length of tag sequence in the ends of primers",
)

# Software group
soft_parser = argparse.ArgumentParser(add_help=False)

soft_group = soft_parser.add_argument_group("software path")

soft_group.add_argument(
    "-vsearch",
    metavar="<STR>",
    help="vsearch path" + "(only needed if vsearch is not in $PATH)",
)

soft_group.add_argument(
    "-threads", metavar="<INT>", default=2, help="threads for vsearch, default=2"
)

soft_group.add_argument(
    "-cid",
    metavar="FLOAT",
    type=float,
    default=0.98,
    dest="cluster_identity",
    help="identity for clustering, default=0.98",
)

## filter group  ##
filter_parser = argparse.ArgumentParser(
    add_help=False,
    description="Use the whole raw" + "dataset (Only adapters should be removed)!",
)

filter_group = filter_parser.add_argument_group("filter arguments")

filter_group.add_argument(
    "-raw",
    metavar="<STR>",
    required=True,
    help="input raw Single-End fastq file, and"
    + " only\nadapters should be removed; supposed on\n"
    + "Phred33 score system (BGISEQ-500)",
)

filter_group.add_argument(
    "-phred",
    metavar="<INT>",
    type=int,
    dest="phred",
    choices=[33, 64],
    default=33,
    help="Phred score system, 33 or 64, default=33"
)

filter_group.add_argument(
    "-e",
    metavar="<INT>",
    type=int,
    dest="expected_err",
    help="expected error threshod, default=10\n"
    + "see more: http://drive5.com/usearch/manual/exp_errs.html",
)

filter_group.add_argument(
    "-q",
    metavar="<INT>",
    type=int,
    dest="quality",
    nargs=2,
    help="filter by base quality; for example: '20 5' means\n"
    + "dropping read which contains more than 5 percent of \n"
    + "quality score < 20 bases.",
)

filter_group.add_argument(
    "-trim",
    dest="trim",
    action="store_true",
    help="whether to trim 5' end of read, it adapts to -e mode\n"
    + "or -q mode",
)
filter_group.add_argument(
    "-n",
    metavar="<INT>",
    type=int,
    default=1,
    help="remove reads containing [INT] Ns, default=1",
)

# ------------------------------------------------------------------------------

## assign group ##
assign_parser = argparse.ArgumentParser(
    add_help=False,
    description="assing clean reads to"
    + "samples by unique tag sequence"
    +"with 100% similarity",
)

assign_group = assign_parser.add_argument_group("assign arguments")

assign_group.add_argument(
    "-primer",
    metavar="<STR>",
    required=True,
    help="taged-primer list, on following format:\n"
    + "Rev001   AAGCTAAACTTCAGGGTGACCAAAAAATCA\n"
    + "For001   AAGCGGTCAACAAATCATAAAGATATTGG\n"
    + "...\n"
    + "this format is necessary!",
)

assign_group.add_argument(
    "-outdir",
    metavar="<STR>",
    default="assigned",
    help="output directory for assignment," + 'default="assigned"',
)

assign_group.add_argument(
    "-tmis",
    metavar="<INT>",
    type=int,
    dest="tag_mismatch",
    default=0,
    help="mismatch number in tag when demultiplexing, default=0",
)

assign_group.add_argument(
    "-pmis",
    metavar="<INT>",
    type=int,
    dest="primer_mismatch",
    default=1,
    help="mismatch number in primer when demultiplexing, default=1",
)

## only assign need
only_assign_parser = argparse.ArgumentParser(add_help=False)
only_assign_group = only_assign_parser.add_argument_group(
    "when only run assign arguments"
)

only_assign_group.add_argument(
    "-fq", metavar="<STR>", required=True, help="cleaned fastq file (*.fq.gz, *.fq)"
)

# ------------------------------------------------------------------------------

## assembly group ##
assembly_parser = argparse.ArgumentParser(
    description="Due to connect HIFI barcode sequence by overlaping two"
    + "consensus sequences which generated from clustering method,"
    + "in this program I use VSEARCH.  You can define the length of overlap,"
    + "how many reads used to make clusters, and whether to check codon"
    + "translation for PCG.",
    add_help=False,
)

assembly_group = assembly_parser.add_argument_group("assembly arguments")

assembly_group.add_argument(
    "-min",
    metavar="INT",
    type=int,
    default=80,
    dest="min_overlap",
    help="minimun length of overlap, default=80",
)

assembly_group.add_argument(
    "-max",
    metavar="INT",
    type=int,
    default=90,
    dest="max_overlap",
    help="maximum length of overlap, default=90",
)

assembly_group.add_argument(
    "-oid",
    metavar="FLOAT",
    type=float,
    default=0.95,
    dest="overlap_identity",
    help="minimun similarity of overlap region, default=0.95",
)

assembly_group.add_argument(
    "-tp",
    metavar="INT",
    type=int,
    dest="cluster_number_needKeep",
    help="how many clusters will be used in" + "assembly, recommend 2",
)

assembly_group.add_argument(
    "-ab",
    metavar="INT",
    type=int,
    dest="abundance_threshod",
    help="keep clusters to assembly if its abundance >=INT ",
)

assembly_group.add_argument(
    "-seqs_lim",
    metavar="INT",
    type=int,
    default=0,
    help="reads number limitation. by default,\n" + "no limitation for input reads",
)

assembly_group.add_argument(
    "-len",
    metavar="INT",
    type=int,
    default=400,
    dest="standard_length",
    help="standard read length, default=400",
)

assembly_group.add_argument(
    "-ds",
    dest="drop_short_read",
    action="store_true",
    help="drop short reads away before assembly",
)

assembly_group.add_argument(
    "-mode",
    metavar="INT",
    type=int,
    choices=[1, 2],
    default=1,
    help="1 or 2; modle 1 is to cluster and keep\n"
    + "most [-tp] abundance clusters, or clusters\n"
    + "abundance more than [-ab], and then make a \n"
    + "consensus sequence for each cluster.\n"
    + "modle 2 is directly to make only one consensus\n"
    + "sequence without clustering. default=1",
)

assembly_group.add_argument(
    "-rc",
    dest="reads_check",
    action="store_true",
    help="whether to check amino acid\n"
    +"translation for reads, default not",
)

# translation need
trans_parser = argparse.ArgumentParser(add_help=False)
trans_group = trans_parser.add_argument_group(
    "translation arguments(when set -rc or -cc)"
)

trans_group.add_argument(
    "-codon",
    metavar="INT",
    type=int,
    dest="codon_table",
    default=5,
    help="codon usage table used to check" + "translation, default=5",
)

trans_group.add_argument(
    "-frame",
    metavar="INT",
    type=int,
    choices=[0, 1, 2],
    default=1,
    help="start codon shift for amino acid" + "translation, default=1",
)

## only assembly need
only_assembly_parser = argparse.ArgumentParser(add_help=False)
only_assembly_group = only_assembly_parser.add_argument_group(
    "only run assembly arguments(not all)"
)

only_assembly_group.add_argument(
    "-list",
    metavar="FILE",
    type=str,
    required=True,
    help="input file, fastq file list. [required]",
)

polish_parser = argparse.ArgumentParser(
    description="polish all assemblies, \n"
    + "to make a confident COI barcode"
    + "reference.",
    add_help=False,
)
polish_group = polish_parser.add_argument_group(
    "polish arguments"
)

polish_group.add_argument(
    "-i",
    metavar="STR",
    type=str,
    dest="coi_input",
    required=True,
    help="COI barcode assemblies (fasta)",
)

polish_group.add_argument(
    "-o",
    metavar="STR",
    type=str,
    dest="coi_output",
    help="polished COI barcode assemblies (fasta)",
)

polish_group.add_argument(
    "-tag",
    metavar="STR",
    type=str,
    dest="sampleMark",
    help="add a mark for each sampel, like: >MARK_001;xxx",
)

polish_group.add_argument(
    "-cc",
    dest="coi_check",
    action="store_false",
    help="whether to check final COI contig's\n"
    + "amino acid translation, default yes",
)

polish_group.add_argument(
    "-cov",
    metavar="INT",
    type=int,
    dest="min_coverage",
    default=5,
    help="minimun coverage of 5' or 3' end allowed, default=5",
)

polish_group.add_argument(
    "-l",
    metavar="INT",
    type=int,
    dest="min_length",
    default=711,
    help="minimun length (with tag and primer) of COI barcode allowed, default=711",
)

polish_group.add_argument(
    "-L",
    metavar="INT",
    type=int,
    dest="max_length",
    default=719,
    help="maximun length (with tag and primer) of COI barcode allowed, default=719",
)

# ------------------------------------------------------------------------------------------------

###############################################################################
#####----------------------- main subcommand parsers --------------------######

description = """

Description

    An automatic pipeline for HIFI-SE400 project, including filtering
    raw reads, assigning reads to samples, assembly HIFI barcodes
    (COI sequences), polished assemblies, and do tax identification.
    See more: https://github.com/comery/HIFI-barcode-SE400

Versions

    1.0.5 (20190409)

Authors

    yangchentao at genomics.cn, BGI.
    mengguanliang at genomics.cn, BGI.
"""

parser = argparse.ArgumentParser(
    prog="HIFI-SE",
    description=description,
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument(
    "-v", "--version",
    action="version",
    version="%(prog)s 1.0.5"
)

subparsers = parser.add_subparsers(dest="command")

########## subcommmands ###########

## all subcommand
parser_all = subparsers.add_parser(
    "all",
    parents=[
        common_parser,
        index_parser,
        soft_parser,
        filter_parser,
        assign_parser,
        assembly_parser,
        trans_parser,
    ],
    formatter_class=argparse.RawTextHelpFormatter,
    help="run filter, assign and assembly.",
)

## filter subcommand
parser_filter = subparsers.add_parser(
    "filter",
    parents=[common_parser, filter_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="remove or trim reads with low quality.",
)

## assign subcommand
parser_assign = subparsers.add_parser(
    "assign",
    parents=[common_parser,
             index_parser,
             only_assign_parser,
             assign_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="assign reads to samples by tags.",
)

## assembly subcommand
parser_assembly = subparsers.add_parser(
    "assembly",
    parents=[common_parser,
             index_parser,
             only_assembly_parser,
             soft_parser,
             assembly_parser,
             trans_parser],
    formatter_class=argparse.RawTextHelpFormatter,
    help="do assembly from assigned reads,\noutput raw HIFI barcodes.",
)

## polish subcommand
parser_polish = subparsers.add_parser(
    "polish",
    parents=[polish_parser,
            index_parser,
            trans_parser,],
    formatter_class=argparse.RawTextHelpFormatter,
    help="polish COI barcode assemblies,\n"
    + "output confident barcodes."
)

## BOLD_identification
parser_bold = subparsers.add_parser(
    "taxonomy",
    parents=[],
    formatter_class=argparse.RawTextHelpFormatter,
    help="do taxa identification on BOLD system\n",
)

###############################################################################
#####---------------------- program execution start ----------------------#####

# -----------------------BOLD identification----------------------#
if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

if sys.argv[1] == "taxonomy":
    # if args.command == 'bold_identification':
    sys.argv = sys.argv[1:]
    sys.exit(bold_identification())

args = parser.parse_args()

# -----------------------arguments checking-----------------------#
## softwares and databases
def check_program_involed(cmd):
    '''
    check program involed whether is executable!
    '''
    result = (
        subprocess.call(
            "type %s" % cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        == 0
    )
    if result:
        return 0
    else:
        print(cmd + " not found!", file=sys.stderr)
        return 1


def files_exist_0_or_1(filelist):
    '''
    check files involed whether are existing!
    '''
    NUM = 0
    for file in filelist:
        if os.path.exists(file):
            NUM += 1
        else:
            print("%s doesn't exist!" % file, file=sys.stderr)
    if len(filelist) == NUM:
        return 0
    else:
        return 1

def print_time(str):
    print(str + " " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
# ----------------------------------------------------------------
## file existing check
errors_found = 0
if args.command == "all":
    errors_found += files_exist_0_or_1([args.raw, args.primer])
elif args.command == "filter":
    errors_found += files_exist_0_or_1([args.raw])
elif args.command == "assign":
    errors_found += files_exist_0_or_1([args.primer])
elif args.command == "assembly":
    errors_found += files_exist_0_or_1([args.list])
elif args.command == "polish":
    errors_found += files_exist_0_or_1([args.coi_input])
else:
    parser.print_help()
    parser.exit()

if args.command in ["all", "assembly"]:
    vsearch = "vsearch"
    if hasattr(args, "vsearch"):
        if args.vsearch:
            vsearch = args.vsearch
    errors_found += check_program_involed(vsearch)

if errors_found > 0:
    parser.exit("Errors found! Exit!")

if hasattr(args, "outpre") and args.outpre.endswith("/"):
    print("outpre is in bad format! no \"/\"")
    exit()

def check_and_open_outhandle(file):
    if os.path.exists(file):
        print("WARRNING: " + file + " exists! now overwriting ...")
    else:
        print("[INFO]: " + "open file " + file + "...")
    out = open(file, 'w')
    return out

# -----------------------functions for filtering------------------#

def parse_se_fastq(fq_fh):
    """
    make generator to read fastq file.
    """
    while True:
        name = fq_fh.readline().strip()
        if len(name) == 0:
            break
        read = fq_fh.readline().strip()
        nothing = fq_fh.readline().strip()
        qual = fq_fh.readline().strip()
        qlist = list(qual)
        yield name, read, qlist

def exp_e(qlist, phred):
    '''
    expected error number(E*) = sum(P), where P is
    the probability that the base call is incorrect.
    '''
    exp = 0
    ascill = [ord(n) - phred for n in qlist]

    for i in ascill:
        exp += 10 ** (-i / 10)

    return exp

def exp_e_trim(seq, qual, E, phred):
    '''
    expected error number(E*) = sum(P), where P is
    the probability that the base call is incorrect.
    the difference between above exp_e funciton is
    that this function will trim the bad based with
    low quality in the 5' end.
    '''
    exp = 0
    good_seq = ''
    good_qual = ''
    for i in range(len(seq)):
        if exp < E:
            good_seq += seq[i]
            good_qual += qual[i]
            ascill_i = ord(qual[i]) - phred
            exp += 10 ** (-ascill_i / 10)
        else:
            return (good_seq, good_qual)
            break
    return (good_seq, good_qual)

def lowquality_rate(qlist, cut_off, phred):
    '''
    calculate the rate of low quality base in read.
    '''
    low_base = 0
    ascill = [ord(n) - phred for n in qlist]

    for i in ascill:
        if i < cut_off:
            low_base += 1

    low_rate = low_base / len(qlist)
    return low_rate

def lowquality_rate_trim(seq, qual, quality_demand, low_rate, phred):
    '''
    calculate the rate of low quality base in read.
    '''
    low_base = 0
    good_seq = ''
    good_qual = ''
    llen = len(seq)
    for i in range(llen):
        ascill_i = ord(qual[i]) - phred
        if low_base < llen * low_rate:
            good_seq += seq[i]
            good_qual += qual[i]
            if ascill_i >= quality_demand:
                # it's a good base
                1
            else:
                low_base += 1
        else:
            return (good_seq, good_qual)
            break

    return (good_seq, good_qual)


# ----------------------functions for assigning-------------------#
def complementation(sequence):
    # make a sequence complement #
    # replace function of string is too low!
    sequence = sequence.upper() ## [a bug fixed], reported by Wu Ping 20181129
    transtable = str.maketrans('ATCG-', 'TAGC-')
    sequence = sequence.translate(transtable)
    return sequence

def comp_rev(sequence):
    # make a sequence complement and reversed #
    sequence = complementation(sequence)
    return sequence[::-1]

def detect_mis(f, r, dict):
    tag_mis = 0
    primer_mis = 0
    strs = [f, r]
    while(strs):
        s1 = strs.pop()
        for s2 in dict.keys():
            if len(s1) == len(s2):
                tag_mis = 0
                primer_mis = 0
                for base in range(len(s1)):
                    if s1[base] is not s2[base]:
                        if base < args.index:
                            tag_mis += 1
                        else:
                            primer_mis += 1
                if (tag_mis <= args.tag_mismatch
                    and primer_mis <= args.primer_mismatch):
                    goal = dict[s2]
                    break
                else:
                    goal = ''
    if len(goal) > 0:
        return goal
    else:
        return False

def dis_barcode(barcode_list):
    # ----count matched bases of two barcodes----#
    dis = []
    while(barcode_list):
        b1 = barcode_list.pop()
        for b2 in barcode_list:
            mismatch = 0
            for base in range(len(b1)):
                if b1[base] is not b2[base]:
                    mismatch += 1
            dis.append(mismatch)
    min_dis = min(dis)
    max_dis = max(dis)
    return (min_dis, max_dis)

# ----------------------functions for assembling------------------#

def comp_rev_list(reads_list):
    # make a list of sequences reverse and complement #
    new_reads_list = []
    for read in reads_list:
        read = comp_rev(read)
        new_reads_list.append(read)

    return new_reads_list

def match(str1, str2):
    # ----count matched bases of two sequences----#
    matched = 0
    for base in range(len(str1)):
        if str1[base] == str2[base]:
            matched += 1
    identity = matched / len(str1)
    return identity

def translate_dnaseq(seq, codon):
    # ---------translate_dnaseq------------#
    l_dna = len(seq)
    if l_dna % 3 is not 0:
        seq = seq[: -(l_dna % 3)]
        # print("your sequence lenght is not tripple" + \
        # "but no worries, I have trimmed well format")
    coding_dna = Seq(seq, generic_dna)
    protein = coding_dna.translate(table=codon)
    if "*" in protein:
        return False
    else:
        return True

def read_fastq(fastq_file, ori):
    # ----------read fastq file-----------#
    if fastq_file.endswith(".gz"):
        fh_file = gzip.open(fastq_file, "rt")
    else:
        fh_file = open(fastq_file, "r")

    reads_count = 0
    good = 0
    short_reads = 0
    seq_checked = []
    # fastq_err = "It's not a correct fastq format.\n"
    for i in parse_se_fastq(fh_file):
        head, sequence, qual = i
        reads_count += 1
        seq_len = len(sequence)
        if args.drop_short_read and seq_len < args.standard_length:
            short_reads += 1
            continue

        # check translation
        if args.reads_check:
            tmp = sequence
            if ori == "f":
                # forward primer length is 25
                tmp_remove = args.index + 25 + args.frame
                tmp = tmp[tmp_remove:]
                if translate_dnaseq(tmp, args.codon_table):
                    seq_checked.append(sequence)
                    good += 1
            else:
                # reverse primer length is 26
                # triming; complementation; triming end; reverse
                tmp_remove = args.index + 26
                tmp = tmp[tmp_remove - 1 :]
                # complementation
                tmp = complementation(tmp)
                # triming end
                if seq_len % 3 is not 0:
                    tmp = tmp[0 : -(seq_len % 3)]
                # reverse #
                tmp = tmp[::-1]
                if translate_dnaseq(tmp, args.codon_table):
                    seq_checked.append(sequence)
                    good += 1
        else:
            # do not check translation #
            seq_checked.append(sequence)
        if args.seqs_lim and reads_count >= args.seqs_lim:
            break

    fh_file.close()
    fh_log.write(ori + ": total reads input: {};".format(reads_count))
    if args.drop_short_read == True:
        fh_log.write(" short reads: {};".format(short_reads))
    if args.reads_check:
        fh_log.write(" check codon good: {}".format(good))
    fh_log.write("\n")

    return seq_checked

def coi_check(contig, codon):
    # ---------------coi_check------------#
    for_trim = args.index + 25 + 1
    rev_trim = args.index + 26
    contig = contig[for_trim:]
    # contig = contig[0:for_trim] #bug!
    contig = contig[:-rev_trim]
    if translate_dnaseq(contig, codon):
        return True
    else:
        return False

def max_length(seqs):
    lens_seqs = [len(n) for n in seqs]
    longest = max(lens_seqs)
    return longest

def depth_table(seqs):
    # ---------------depth_table----------#
    longest = max_length(seqs)
    consensus = ""
    for m in range(longest):
        depth = {}
        for align in seqs:
            if (m + 1) > len(align):
                continue # when sequence is too short, it will be out of index.
            if align[m] in depth.keys():
                depth[align[m]] += 1
            else:
                depth[align[m]] = 1

        if "-" in depth.keys():
            del depth["-"]

        # sorted_base = sorted(depth.iteritems(), key=lambda x: x[1], reverse=True)#
        # sorted_base = sorted(depth, key=depth.__getitem__,reverse=True)
        sorted_base = sorted(depth, key=lambda k: (depth[k], k), reverse=True)

        # for some case, different bases have same deepth #

        consensus += sorted_base[0]

        depth = {}
    return consensus

def report_depth(table, title, seq, read_len, step, ori):
    # report the depth of each site
    depth_sum = {}
    reports = []
    four_bases = ("A", "T", "C", "G")
    reports = table[seq]
    if ori == "f":
        for_rest_len = read_len - step
        for x in range(for_rest_len):
            real_position = x + 1
            fh_depth.write(title + "\t" + str(real_position) + "\t")

            for y in range(len(reports)):
                if reports[y][x] in depth_sum.keys():
                    depth_sum[reports[y][x]] += 1
                else:
                    depth_sum[reports[y][x]] = 1

            for base in four_bases:
                if base in depth_sum.keys():
                    fh_depth.write(base + ":" + str(depth_sum[base]) + "\t")
                else:
                    fh_depth.write(base + ":0" + "\t")
            fh_depth.write("\n")
            depth_sum = {}

    else:
        for x in range(step, read_len):
            real_position = read_len - step + x + 1
            fh_depth.write(title + "\t" + str(real_position) + "\t")

            for y in range(len(reports)):
                if reports[y][x] in depth_sum.keys():
                    depth_sum[reports[y][x]] += 1
                else:
                    depth_sum[reports[y][x]] = 1

            for base in four_bases:
                if base in depth_sum.keys():
                    fh_depth.write(base + ":" + str(depth_sum[base]) + "\t")
                else:
                    fh_depth.write(base + ":0" + "\t")
            fh_depth.write("\n")
            depth_sum = {}


def mode_identical(seqs):
    """
    in mode 1, using all codon checked reads to cluster with 100% identity,
    this is easy to achieve.
    """
    abu = {}
    for s in seqs:
        if s in abu.keys():
            abu[s] += 1
        else:
            abu[s] = 1
    seqs = []

    # sort and cluster with 100% identity#
    # sorted_seq = sorted(abu, key=abu.__getitem__,reverse=True)
    sorted_seq = sorted(abu, key=lambda k: (abu[k], k), reverse=True)

    most_abuns = {}
    # keep top 1#
    if abu[sorted_seq[0]] > 1:
        # keep top 2, if top1 has same abundance with top2's#
        if abu[sorted_seq[0]] == abu[sorted_seq[1]]:
            most_abuns[sorted_seq[0]] = abu[sorted_seq[0]]
            most_abuns[sorted_seq[1]] = abu[sorted_seq[1]]
            print("Alert: top1 and top2 have same number of sequences!\n")
            fh_log.write("Alert: top1 and top2 have same number of sequences!\n")
        else:
            most_abuns[sorted_seq[0]] = abu[sorted_seq[0]]
    else:
        # all reads are unique, so that's wrong! #
        most_abuns[sorted_seq[0]] = abu[sorted_seq[0]]
        print("Alert: all reads in file are uniqe! clustering is meaningless!\n")
    abu = {}

    return most_abuns

def merge_matrix(matrix1, matrix2):
    # merge two matrix
    for i in range(len(matrix2)):
        matrix1.append(matrix2[i])

    return matrix1

def repair_short_reads(reads):
    repaired = []
    lmax = max_length(reads)
    for i in reads:
        if len(i) < lmax:
            # completion with - in the head for end reads
            a = i + "-" * (lmax - len(i))
        else:
            a = i
        repaired.append(a)
    return repaired

def mode_vsearch(seqs):
    """
    in mode 1, but clustering identity is not 100%,
    it is hard to archieve by perl, so I use VSEARCH to make it.
    I just keep three sequences with the highest abundance
    after sorting all clusters.
    HASH table %cluster_tables contains pairs,
    which key = one cluster's consensus sequence, val = all sequences in
    this cluster with 2D array format.
    'consensus sequence' -> 'bases matrix build from same cluster'
    """
    cluster_tables = {}
    pid = os.getpid()
    temp_fasta = "temp.fa" + "." + str(pid)
    temp_uc = "temp.uc" + "." + str(pid)
    with open(temp_fasta, "w") as TM:
        for i in range(len(seqs)):
            TM.write(">" + str(i) + "\n" + seqs[i] + "\n")

    vsearch_cmd = (
        vsearch
        + " --cluster_fast "
        + temp_fasta
        + " --threads "
        + str(args.threads)
        + " --quiet "
        + " --uc "
        + temp_uc
        + " --id "
        + str(args.cluster_identity)
    )
    # print("now run: " + vsearch_cmd)
    subprocess.call(vsearch_cmd, shell=True)
    # to store clusters, clusters[represent ID] = [clustered IDs]
    clusters = {}
    # to store each cluster's abundance, for sorting clusters
    # e.g. count[represent ID] = 100
    count = {}

    # open temp.uc and statistic each cluster's abundance.
    with open(temp_uc, "r") as uc:
        # there are "H","S","C" in the head of line
        for line in uc.readlines():
            if line[0] is not "H":
                continue
            array = line.split()
            if array[9] in clusters.keys():
                clusters[array[9]].append(array[8])
                count[array[9]] += 1
            else:
                clusters[array[9]] = []
                clusters[array[9]].append(array[9])
                count[array[9]] = 1

    # sorting clusters by abundance.#
    # sorted_clusters = sorted(count, key=count.__getitem__,reverse=True)
    sorted_clusters = sorted(count, key=lambda k: (count[k], k), reverse=True)

    if args.cluster_number_needKeep:
        # if set "-tp", keep top N clusters to assembly
        keep = args.cluster_number_needKeep
        sorted_clusters = sorted_clusters[0:keep]
    elif args.abundance_threshod:
        # if set "-ab", keep all clusters of abundance > ab
        while sorted_clusters:
            item = sorted_clusters.pop()
            if count[item] < args.abundance_threshod:
                pass
            else:
                sorted_clusters.append(item)
                break
    elif len(sorted_clusters) > 1:
        # if set nothing, I will set it to -tp 2
        # if second most abundant sequence less than 1/10 of first,
        # remove it!# of course it is just for when tp==2
        sorted_clusters = sorted_clusters[0:2]
        if count[sorted_clusters[1]] < count[sorted_clusters[0]] / 10:
            sorted_clusters.pop()

    for k in sorted_clusters:
        k_seqs = []
        matrix = []
        for s in clusters[k]:
            s = int(s)
            k_seqs.append(seqs[s])
            # split reads into bases-arrays#
        # if your reads are trimed, so they are uneven, this is to
        # normalize length
        k_seqs = repair_short_reads(k_seqs)
        con = depth_table(k_seqs)
        for s in k_seqs:
            matrix.append(list(s))
        '''
        #!!!---BUG---report
        # if make a consensus sequence for each cluster,
        # you may face this problem: different clusters have same
        # consensus sequence, so previous key will be masked.
        # On! fuck!
        #
        # Take it easy! I have a solution!
        # why not to merge two pairs which have same keys.
        '''
        if con in cluster_tables.keys():
            # merge this matric and previous same matrix#
            flash = merge_matrix(cluster_tables[con], matrix)
            cluster_tables[con] = flash
        else:
            cluster_tables[con] = matrix

    return cluster_tables

def mode_consensus(seqs):
    # In mode 2, this is only table.#
    matrix = []
    cons_table = {}
    for s in seqs:
        matrix.append(list(s))

    consensus = depth_table(seqs)
    cons_table[consensus] = matrix
    return cons_table

def addtwodimdict(thedict, key_a, key_b, val):
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})

# ------------------------filter process--------------------------#
if args.command in ["all", "filter"]:

    print_time("[INFO]: Filtering start:")

    filtered_outfile = args.outpre + "_filter_highqual.fastq"
    out = check_and_open_outhandle(filtered_outfile)

    logfile = args.outpre + "_filter_log.txt"
    log = check_and_open_outhandle(logfile)

    # ini state
    total = 0
    clean = 0
    nn = 0
    if args.trim == True:
        read_lens = []
        min_read_len = 0
        max_read_len = 0
        trimed_base = 0
    else:
        err = open(args.outpre + "_filter_lowqual.fastq", "w")

    if args.expected_err and args.quality:
        print(
            "[Bad Arguments]:"
            + " -e argument is confilicting with -q,"
            + " can not using in the same time"
        )
        exit()
    elif args.quality:
        high_qual = args.quality[0]
        low_qual_cont = args.quality[1] / 100
        FILTER_TYPE = 2
        log.write("Filtering by quality score: {}, content: {}%".format(
            args.quality, args.quality[1])
            + "\n")
    else:
        FILTER_TYPE = 1
        if not args.expected_err:
            args.expected_err = 10
        log.write("Filtering by expected_err: {}".format(args.expected_err)
                  + "\n")

    if args.raw.endswith(".gz"):
        fh = gzip.open(args.raw, "rt")
    else:
        fh = open(args.raw, "r")

    for i in parse_se_fastq(fh):
        head, seq, qual = i
        qual_str = "".join(qual)
        total += 1
        N_count = seq.count("N")

        if N_count < args.n:
            if FILTER_TYPE == 1:
                if args.trim == True:
                    (good_seq, good_qual) = exp_e_trim(seq,
                                                       qual,
                                                       args.expected_err,
                                                       args.phred)
                    out.write(head
                              + "\n"
                              + good_seq
                              + "\n+\n"
                              + good_qual
                              + "\n")
                    read_lens.append(len(good_seq))
                    clean += 1
                elif exp_e(qual, args.phred) <= args.expected_err:
                    out.write(head + "\n" + seq + "\n" + "+\n" + qual_str + "\n")
                    clean += 1
                else:
                    err.write(head + "\n" + seq + "\n" + "+\n" + qual_str + "\n")
            else:
                # filter_type == 2
                if args.trim == True:
                    (good_seq, good_qual) = lowquality_rate_trim(seq, qual,
                                                                 high_qual,
                                                                 low_qual_cont,
                                                                args.phred)
                    out.write(head
                              + "\n"
                              + good_seq
                              + "\n+\n"
                              + good_qual
                              + "\n")
                    read_lens.append(len(good_seq))
                    clean += 1
                elif lowquality_rate(qual, high_qual, args.phred) > low_qual_cont:
                    out.write(head
                              + "\n"
                              + seq
                              + "\n+\n"
                              + qual_str
                              + "\n")
                    clean += 1
                else:
                    err.write(head
                              + "\n"
                              + seq
                              + "\n+\n"
                              + qual_str
                              + "\n")
        else:
            nn += 1
    fh.close()

    log.write("total reads:\t{}".format(total) + "\n")
    log.write("clean reads:\t{}".format(clean) + "\n")
    log.write("containing N reads:\t{}".format(nn) + "\n")
    if args.trim == True:
        for le in read_lens:
            trimed_base += le
        average_trimed = trimed_base / clean
        min_read_len = min(read_lens)
        max_read_len = max(read_lens)
        log.write("clean base:\t{}".format(trimed_base) + "\n")
        log.write("min read length:\t{}".format(min_read_len) + "\n")
        log.write("max read length:\t{}".format(max_read_len) + "\n")
        log.write("average of trimed read:\t{0:.2f}".format(average_trimed))
    else:
        err.close()

    log.close()
    out.close()

    print_time("[INFO]: Filtering done:")

# ------------------------assign process--------------------------#

if args.command in ["all", "assign"]:
    print_time("[INFO]: Assigning start:")

    if args.command == "all":
        args.fq = filtered_outfile
        assigned_outdir = os.path.abspath(args.outpre + "_assign")
        ErrFile = open(assigned_outdir + "_assign_err.fasta", "w")

    elif args.command == "assign":
        assigned_outdir = os.path.abspath(args.outdir)
        ErrFile = open(assigned_outdir + "_err.fasta", "w")

    indexlen = args.index

    if os.path.exists(assigned_outdir) == False:
        os.mkdir(assigned_outdir)

    pris = {}
    indp = {}
    FH = {}
    barcodes = []

    with open(args.primer, "r") as p:
        primer_lines = 0
        for i in p.readlines():
            primer_lines += 1
            i = i.strip()
            arr = i.split()
            if len(arr) != 2:
                print("[ERROR]: primer set is not well-formated")
                exit()
            sam = arr[0]
            ipr = arr[1]
            FH[ipr] = arr[0]
            if sam in indp.keys():
                print("[INFO]: " + arr[0] + "show twice in primer set")
            else:
                ori = sam[0:3]
                num = sam[-3:]

                pris[num] = {}
                pris[num][ori] = ipr
                indp[sam] = ipr

            if "For" in sam:
                plenf = len(ipr)
                primerF = ipr[indexlen:]
                # save barcodes in a array, analyze later
                barcodes.append(ipr[0:indexlen])

            if "Rev" in sam:
                plenr = len(ipr)
                primerR = ipr[indexlen:]
    # check primer lines
    if primer_lines % 2 != 0:
        print(
            "[ERROR]: primer lines ({}) is wrong!".format(primer_lines)
        )
        print(
            "[INFO]: the primer file need to have each"
            + "forward and reverse primer"
        )
        exit()

    # analysis barcodes and demultiplex argument(mismatch)
    (min_dis, max_dis) = dis_barcode(barcodes)
    print("[INFO]: min distance among barcodes is {}".format(min_dis))
    print("[INFO]: max distance among barcodes is {}".format(max_dis))
    if args.tag_mismatch and args.tag_mismatch > (min_dis - 1):
        print("[ERROR]: mismatch you set is too large to demultiplex,"
              + " it must be smaller than min distance ({})".format(min_dis))
        exit()
    # ini dict of count_total and count_assigned for statistic.
    count_assigned = {}
    for fh in FH.values():
        count_assigned[fh] = 0


    neg_priF = comp_rev(primerF)
    neg_priR = comp_rev(primerR)
    assigned_list =  args.outpre + "_assign.list"

    with open(assigned_list, "w") as ls:
        sorted_sample = sorted(pris.keys())
        for s in sorted_sample:
            ls.write(
                assigned_outdir
                + "/For"
                + s
                + ".fastq"
                + "\t"
                + assigned_outdir
                + "/Rev"
                + s
                + ".fastq"
                + "\n"
            )
    # open all assigned files
    filehandle = {}
    for sam in indp.keys():
        filehandle[sam] = open(assigned_outdir + "/" + sam + ".fastq", "w")

    seqnum = 0
    err = 0
    assigned = 0

    if args.fq.endswith(".gz"):
        fh = gzip.open(args.fq, "rt")
    else:
        fh = open(args.fq, "r")
    for i in parse_se_fastq(fh):
        head, seq, qual = i
        qual_str = "".join(qual)
        seqnum += 1
        headf = seq[0:plenf]
        headr = seq[0:plenr]
        # cut head (max of for and rev) to make a tmp sequence
        len_head_cut = max(plenf, plenr)
        tmp = seq[len_head_cut:]
        # if primer in the wrong position, remove this reads.
        if (primerF in tmp or primerR in tmp
            or neg_priF in tmp or neg_priR in tmp):
            ErrFile.write(">" + str(seqnum) + "_priErr\n" + seq + "\n")
            err += 1

        elif headf in FH.keys():
            count_assigned[FH[headf]] += 1
            filehandle[FH[headf]].write(
                "@" + FH[headf] + "_" + str(seqnum) + "\n" + seq + "\n"
                + "+\n" + qual_str + "\n")
            assigned += 1

        elif headr in FH.keys():
            count_assigned[FH[headr]] += 1
            filehandle[FH[headr]].write(
                "@" + FH[headr] + "_" + str(seqnum) + "\n" + seq + "\n"
                + "+\n" + qual_str + "\n")
            assigned += 1
        else:
            # much more likely to be a mismatch
            potential_target = detect_mis(headf, headr, FH)
            if potential_target:
                filehandle[potential_target].write(
                    "@"
                    + potential_target
                    + "_"
                    + str(seqnum)
                    + "\n"
                    + seq
                    + "\n+\n"
                    + qual_str
                    + "\n")
                assigned += 1
                count_assigned[potential_target] += 1
            else:
                err += 1
                ErrFile.write(">" + str(seqnum) + "\n" + seq + "\n")
    #close args.fq
    fh.close

    ErrFile.close()
    # close all assigned files
    for fh in filehandle.values():
        fh.close()

    # report assignment information
    with open(args.outpre + "_assign.log", "w") as log:
        log.write("total reads:\t{}\n".format(seqnum))
        log.write("err reads:\t{}\n".format(err))
        log.write("assigned:\t{}\n".format(assigned))
        for i in sorted(count_assigned.keys()):
            log.write(
                i + "\t" + str(count_assigned[i]) + "\n"
            )

    print_time("[INFO]: Assigning done:")

# ------------------------assembly process--------------------------#
if args.command in ["all", "assembly"]:
    print_time("[INFO]: Assembling start:")

    if args.cluster_number_needKeep and args.abundance_threshod:
        print(
            "Bad arguments:\n\t"
            + "-tp argument is confilicting with -ab,"
            + " can not using in the same time"
        )
        exit()

    if args.min_overlap > 83:
        print("[ERROR]: "
            + "For COI barcodes, by and large, overlaped length is 83 bp, so\
              {} is not proper!".format(args.min_overlap)
        )
        exit()
    if args.max_overlap < args.min_overlap:
        print("[ERROR]: maximum overlap length must be large than minimun")
        exit()

    if args.command == "all":
        args.list = assigned_list  # list generated from assign step

    assembly_result = args.outpre + "_assembly.fasta"
    fh_out = check_and_open_outhandle(assembly_result)
    fh_log = check_and_open_outhandle(args.outpre + "_assembly.log")
    fh_depth = check_and_open_outhandle(args.outpre + "_assembly.depth")

    fh_log.write("## assigned reads list file = " + args.list + "\n")

    if args.seqs_lim:
        fh_log.write("## reads input limitation: " + str(args.seqs_lim) + "\n")
    else:
        fh_log.write("## Using all reads to make consensus, no limitation\n")

    fh_log.write("## input read length = " + str(args.standard_length) + "\n")
    fh_log.write("## vsearch path = " + vsearch + "\n")

    if args.reads_check:
        fh_log.write("## check codon translation = yes\n")
    else:
        fh_log.write("## check codon translation = no\n")

    fh_log.write("## consensus mode = " + str(args.mode) + "\n")

    if args.mode == 1:
        fh_log.write("## clustering identity = "
                     + str(args.cluster_identity)
                     + "\n")

    fh_log.write("## overlaping identity = "
                 + str(args.overlap_identity)
                 + "\n")
    fh_log.write("## min overlap = " + str(args.min_overlap) + "\n")
    fh_log.write("## max overlap = " + str(args.max_overlap) + "\n")

    list_format_info = (
        "[ERROR]: your list is not well formated!"
        + "for example:\n\t"
        + "/path/test_For001.fastq\t/path/test_Rev001.fastq"
    )

    # --------------main-----------------------#
    barcodes_count = 0
    run_again = []
    run_again_file = ""
    try:
        with open(args.list) as fh_list:
            lines = fh_list.readlines()
    except FileNotFoundError:
        print("[ERROR]: can not find " + args.list)
        exit(0)

    for line in lines:
        success_or_not = False
        line = line.rstrip()
        tmp_list = line.split()
        if len(tmp_list) != 2:
            print(list_format_info)
            exit(0)
        else:
            forward = tmp_list[0]
            reverse = tmp_list[1]

            for_name = os.path.basename(forward).split(".")[0]
            rev_name = os.path.basename(reverse).split(".")[0]
            outname = for_name + "_" + rev_name
            short_outname = outname[-3:]
            fh_log.write("//processing " + outname + "\n")
            # if file is empty, continue
            if os.path.getsize(forward) == 0 or os.path.getsize(reverse) == 0:
                fh_log.write("! Eithor Forward or Reverse file is empty!")
                fh_log.write("\n")
                continue

            seq_checked_for = read_fastq(forward, "f")
            seq_checked_rev = read_fastq(reverse, "r")

            if len(seq_checked_for) == 0 or len(seq_checked_rev) == 0:
                fh_log.write("! Eithor Forward or Reverse reads is NONE!")
                fh_log.write("\n")
                run_again.append(short_outname)
                run_again_file += line + "\n"
                continue

            # here table_f and table_r are two quotes of two dicts,
            # key is consensus sequence value is a 2D array,
            # including each cluster's bases.

            if args.mode == 1 and args.cluster_identity == 1:
                table_f = mode_identical(seq_checked_for)
                seq_checked_rev = comp_rev_list(seq_checked_rev)
                table_r = mode_identical(seq_checked_rev)

            elif args.mode == 1 and args.cluster_identity < 1:
                table_f = mode_vsearch(seq_checked_for)
                seq_checked_rev = comp_rev_list(seq_checked_rev)
                table_r = mode_vsearch(seq_checked_rev)

            else:
                # mod == 2
                # in mod2, though there is only one sequence, also using a array to store.

                seq_checked_for = repair_short_reads(seq_checked_for)
                table_f = mode_consensus(seq_checked_for)
                seq_checked_rev = repair_short_reads(seq_checked_rev)
                seq_checked_rev = comp_rev_list(seq_checked_rev)
                table_r = mode_consensus(seq_checked_rev)

            # sort array by it's abundance
            consensus_for = sorted(
                table_f.keys(), key=lambda x: len(table_f[x]), reverse=True
            )
            consensus_rev = sorted(
                table_r.keys(), key=lambda x: len(table_r[x]), reverse=True
            )

            ##-----------anchoring overlap site--------#
            i = 0
            j = 0
            len_conF = len(consensus_for)
            len_conR = len(consensus_rev)
            for i in range(0, len_conF):
                cluster_f = consensus_for[i]
                pos1 = i + 1
                abundance_f = len(table_f[cluster_f])

                fh_log.write(
                    ">"
                    + short_outname
                    + "_f_"
                    + str(pos1)
                    + "\tsize="
                    + str(abundance_f)
                    + "\n"
                    + cluster_f
                    + "\n"
                )

                for j in range(0, len_conR):
                    cluster_r = consensus_rev[j]
                    pos2 = j + 1
                    abundance_r = len(table_r[cluster_r])

                    # print just once
                    if pos1 == 1:
                        fh_log.write(
                            ">"
                            + short_outname
                            + "_r_"
                            + str(pos2)
                            + "\tsize="
                            + str(abundance_r)
                            + "\n"
                            + cluster_r
                            + "\n"
                        )
                    c_size = (abundance_f + abundance_r) / 2
                    c_size = str(c_size)

                    read0 = cluster_f[-args.max_overlap :]
                    read1 = cluster_r[0 : args.max_overlap]

                    singal = 0
                    overlaps = {}
                    for s in range(args.min_overlap, args.max_overlap + 1):
                        l0 = read0[-s:]
                        l1 = read1[0:s]
                        tmp_identity = match(l0, l1)
                        if tmp_identity == 1:
                            overlaps[s] = 1
                            fh_log.write(
                                str(pos1)
                                + "-"
                                + str(pos2)
                                + " ---> find position "
                                + str(s)
                                + " with 100% mathch\n"
                            )
                            # find best result, so exit loop #
                            break

                        elif tmp_identity >= args.overlap_identity:
                            overlaps[s] = tmp_identity

                    # find best overlaping result in all potenial positions
                    # candidates = sorted(overlaps.items(),
                    # lambda x, y: cmp(x[1], y[1]), reverse=True)
                    candidates = sorted(
                        overlaps, key=overlaps.__getitem__, reverse=True
                    )

                    if len(candidates) > 0:
                        success_or_not = True
                        potenial = candidates[0]
                        s0 = read0[-potenial:]
                        s1 = read1[0:potenial]

                        if args.mode == 1 and args.cluster_identity == 1:
                            """
                            if cid == 1, so depth of forward and reverse sequence
                            is value of hash(table_f,table_r)
                            otherwise, mod==1,cid<1 or mod==2 have sample hash structure.
                            """
                            fh_depth.write("# cid = 1, no need to report depth!")

                            if table_f[cluster_f] > table_r[cluster_r]:
                                makeup_consensus = (
                                    cluster_f
                                    + cluster_r[potenial - args.standard_length :]
                                )
                            elif table_f[cluster_f] < table_r[cluster_r]:
                                makeup_consensus = (
                                    cluster_f[: args.standard_length - potenial]
                                    + cluster_r
                                )
                            else:
                                makeup_consensus = (
                                    cluster_f
                                    + cluster_r[potenial - args.standard_length :]
                                )
                                fh_log.write(
                                    "forward and reverse have same depth, \
                                    use forward region\n"
                                )

                        else:
                            correct = ""
                            # report forward depth
                            title = short_outname + "_" + str(pos1) + "-" + str(pos2)
                            report_depth(
                                table_f,
                                title,
                                cluster_f,
                                args.standard_length,
                                potenial,
                                "f",
                            )
                            fh_depth.write(
                                short_outname
                                + "_"
                                + str(pos1)
                                + "-"
                                + str(pos2)
                                + "\t\t-----overlap start-----\n"
                            )

                            # compare each base from forward and reverse to keep one
                            # forward == reverse
                            # forward ne reverse, depth(forward) > depth(reverse)
                            # forward ne reverse, depth(forward) < depth(reverse)

                            for p in range(len(s0)):
                                # site is changed, be careful!#
                                tmp_loca0 = args.standard_length - potenial + p
                                # 2D array
                                for_tab = table_f[cluster_f]
                                rev_tab = table_r[cluster_r]
                                (for_depth, rev_depth) = (0, 0)
                                sum = {}
                                sum[s0[p]] = 0
                                sum[s1[p]] = 0

                                for x in range(len(for_tab)):
                                    if for_tab[x][tmp_loca0] == s0[p]:
                                        for_depth += 1
                                        sum[s0[p]] += 1

                                for y in range(len(rev_tab)):
                                    if rev_tab[y][p] == s1[p]:
                                        rev_depth += 1
                                        sum[s1[p]] += 1

                                tmp_loca1 = tmp_loca0 + 1
                                fh_depth.write(
                                    short_outname
                                    + "_"
                                    + str(pos1)
                                    + "-"
                                    + str(pos2)
                                    + "\t"
                                    + str(tmp_loca1)
                                    + "\t"
                                )
                                four_bases = ("A", "T", "C", "G")
                                for base in four_bases:
                                    if base in sum.keys():
                                        fh_depth.write(
                                            base + ":" + str(sum[base]) + "\t"
                                        )
                                    else:
                                        fh_depth.write(base + ":0" + "\t")
                                # empty dict sum for next round #
                                sum = {}

                                if s0[p] == s1[p]:
                                    correct += s0[p]
                                    info = s0[p] + "=" + s1[p]
                                else:
                                    if for_depth > rev_depth:
                                        correct += s0[p]
                                        info = (
                                            "[F="
                                            + s0[p]
                                            + ":"
                                            + str(for_depth)
                                            + " > R="
                                            + s1[p]
                                            + ":"
                                            + str(rev_depth)
                                            + "]"
                                        )
                                    elif for_depth == rev_depth:
                                        fh_log.write(
                                            short_outname
                                            + "\t"
                                            + str(tmp_loca1)
                                            + "\tcould be a heterozygosis site!\n"
                                        )
                                        correct += s0[p]
                                        info = (
                                            "[F="
                                            + s0[p]
                                            + ":"
                                            + str(for_depth)
                                            + " = R="
                                            + s1[p]
                                            + ":"
                                            + str(rev_depth)
                                            + "]"
                                        )
                                    else:
                                        correct += s1[p]
                                        info = (
                                            "[F="
                                            + s0[p]
                                            + ":"
                                            + str(for_depth)
                                            + " < R="
                                            + s1[p]
                                            + ":"
                                            + str(rev_depth)
                                            + "]"
                                        )

                                fh_depth.write(info + "\n")

                            fh_depth.write(
                                short_outname
                                + "_"
                                + str(pos1)
                                + "-"
                                + str(pos2)
                                + "\t\t-----overlap end-----\n"
                            )
                            # report reverse depth #
                            title = (
                                short_outname + "_" + str(pos1) + "-" + str(pos2) + "\t"
                            )
                            report_depth(
                                table_r,
                                title,
                                cluster_r,
                                args.standard_length,
                                potenial,
                                "r",
                            )
                            makeup_consensus = (
                                cluster_f[: args.standard_length - potenial]
                                + correct
                                + cluster_r[potenial - args.standard_length :]
                            )

                        len_makeup_consensus = str(len(makeup_consensus))
                        this_oid = overlaps[potenial] * 100
                        this_oid = str("%.2f" % this_oid)

                        # write into output file #
                        barcodes_count += 1
                        fh_out.write(
                            ">"
                            + short_outname
                            + "_"
                            + str(pos1)
                            + "-"
                            + str(pos2)
                            + ";"
                            + str(abundance_f)
                            + "_"
                            + str(abundance_r)
                            + ";size="
                            + c_size
                            + ";overPos="
                            + str(potenial)
                            + ";oid="
                            + this_oid
                            + "%"
                            + ";len="
                            + len_makeup_consensus
                            + "\n"
                            + makeup_consensus
                            + "\n"
                        )
                        # if check result, and ok so write into output_checked file #

                    else:
                        fh_log.write(
                            "<<" + str(pos1) + "-" + str(pos2) + " has no result!\n"
                        )
        if success_or_not == False:
            run_again.append(short_outname)
            run_again_file += line + "\n"
    fh_out.close()

    if args.mode == 1 and args.cluster_identity < 1:
        rm_tmp_cmd = "rm temp.fa.* temp.uc.*"
        os.system(rm_tmp_cmd)

    print("[INFO]: Total barcodes generated: {}".format(barcodes_count))
    print_time("[INFO]: Assembling done:")
    print("\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

    # whether to run second round of assembly
    if len(run_again) > 0:
        rerun_list = args.outpre + "_rerun.list"
        print(
            "Kindly recommend you to run these samples in\n"
            + rerun_list
            + " again, with -rc option,\nor with a larger number of [-tp],"
            + " or try mode 2:"
        )
        print("------")
        with open(args.outpre + "_rerun.list", "w") as rr:
            rr.write(run_again_file)
            for r in run_again:
                print(r)
        print("------")
        print("And, finally run: 'HIFI-SE.py polish' to\npolish assemblies.")
        print(">>>>>>>>>>>>>>>>>>>>END>>>>>>>>>>>>>>>>>>>>>\n")

#--------------------polish process-------------------------------------------
if args.command == "polish":
    from Bio import SeqIO
    if args.coi_output:
        polish_outfile = args.coi_output
    else:
        polish_outfile = args.coi_input + ".polished"

    if args.sampleMark:
        marker = args.sampleMark
    else:
        marker = ""

    coiout = open(polish_outfile,'w')
    with open(args.coi_input, "r") as handle:
        HASH_sam_abu = {}
        ARR_sam_abu = {}
        for record in SeqIO.parse(handle, "fasta"):
            item_id = record.id.split(";")
            sample_tag = item_id[0].split("_")
            sample = sample_tag[0]
            coverages = item_id[1].split("_")
            for_coverage = int(coverages[0])
            rev_coverage = int(coverages[1])
            abu = (for_coverage + rev_coverage) / 2
            length = len(record.seq)
            # remove low record with low coverage or short length
            if (for_coverage < args.min_coverage
                or rev_coverage < args.min_coverage
                or length < args.min_length
                or length > args.max_length):
                continue

            elif args.coi_check and coi_check(str(record.seq), args.codon_table) == False:
                print(str(sample_tag) + " translation failed")
                continue

            else:
                addtwodimdict(HASH_sam_abu, sample, abu, record.seq)
                if sample in ARR_sam_abu.keys():
                    ARR_sam_abu[sample].append(abu)
                else:
                    ARR_sam_abu[sample] = []
                    ARR_sam_abu[sample].append(abu)

    # pick up most confident COI barcode.
    sorted_samples = sorted(ARR_sam_abu.keys())
    polished_count = len(sorted_samples)
    for s in sorted_samples:
        abus = sorted(ARR_sam_abu[s], reverse=True)
        top_abu = abus[0]
        seq = HASH_sam_abu[s][top_abu]
        seqlen = len(seq)
        coiout.write(
            ">"
            + marker
            + "_"
            + str(s)
            + ";size="
            + str(top_abu)
            + ";"
            + "len="
            + str(seqlen)
            + "\n"
            + str(seq)
            + "\n"
        )
    coiout.close()
    print("[INFO]: Total barcodes polished: {}".format(polished_count))

print("[INFO]: total run time: {0:.2f}".format(time.time() - t) + "s")
