#!/usr/bin/env python3
import os
import sys
import time
import gzip
import argparse
import subprocess
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from bold_identification.BOLDv4_identification_selenium import main as bold_identification

t = time.time()

try:
	import Bio
except:
	sys.exit("package biopython not found! Please install it!")


###############################################################################
#####------------------------- parameters --------------------------------#####

## common group  ##
common_parser = argparse.ArgumentParser(add_help=False)

common_group = common_parser.add_argument_group('common arguments')

common_group.add_argument("-outpre", metavar="<STR>", required=True,
						help="outprefix for process.")
## index group ##
index_parser = argparse.ArgumentParser(add_help=False)

index_group = index_parser.add_argument_group('index arguments')

index_group.add_argument('-index', metavar='INT', type=int, required=True,
						help="index sequence lenght")

# Software group
soft_parser = argparse.ArgumentParser(add_help=False)

soft_group = soft_parser.add_argument_group('software path')

soft_group.add_argument("-vsearch", metavar="<STR>", help="vsearch path" +\
						"(only needed if vsearch is not in $PATH)")

soft_group.add_argument("-threads", metavar="<INT>", help="threads for vsearch",
						default=2)

soft_group.add_argument('-cid', metavar='FLOAT', type=float, default=0.98, 
						dest='cluster_identity',
						help="identity for clustering [0.98]")

## filter group  ##
filter_parser = argparse.ArgumentParser(add_help=False,
	description="Use the raw whole dataset (Only adapters should be removed)!")

filter_group = filter_parser.add_argument_group('filter arguments')

filter_group.add_argument("-raw", metavar="<STR>", required=True,
						help="input raw singled-end fastq file, (Phred33)")

filter_group.add_argument("-e", metavar="<INT>", type=int, dest="expected_err",
						help="expected error number threshod, P = 10–Q/10, default=10")

filter_group.add_argument("-q", metavar="<INT>", type=int, dest="quality",nargs=2,
						help="filter by quality method,  Q = –10 log10(P),\n" +\
						"filter out low quality reads. example: 20 5, it means\n" +\
						"dropping read which contains more than 5 percent of \n" +\
						"quality score < 20 bases.")

filter_group.add_argument("-n", metavar="<INT>", type=int, default=1,
						help="remove reads containing [INT] Ns, default=1")

#------------------------------------------------------------------------------

## assign group ##
assign_parser = argparse.ArgumentParser(add_help=False,
	description="assing clean reads to samples by unique tag sequence with 100% similarity")

assign_group = assign_parser.add_argument_group('assign arguments')

assign_group.add_argument("-primer", metavar="<STR>", required=True,
							help="taged primer list, like following lines:\n" +\
							"Rev001	AAGCTAAACTTCAGGGTGACCAAAAAATCA\n" +\
							"For001	AAGCGGTCAACAAATCATAAAGATATTGG\n" +\
							"...\n" +\
							"this format is necessary!")

assign_group.add_argument("-outdir", metavar="<STR>", default='assigned', 
							help="output directory for assignment")
## only assign need
only_assign_parser = argparse.ArgumentParser(add_help=False)
only_assign_group = only_assign_parser.add_argument_group('when only run assign arguments')

only_assign_group.add_argument("-fq", metavar="<STR>", required=True,
							help="cleaned fastq file ")

#------------------------------------------------------------------------------

## assembly group ##
assembly_parser = argparse.ArgumentParser(
	description="Due to connect HIFI barcode sequence by overlaping two " +\
	"consensus sequences which generated from clustering method, " +\
	"in this program I use VSEARCH.  You can define the length of overlap," +\
	"how many reads used to make clusters, and whether to check codon translation for PCG.",
	add_help=False)

assembly_group = assembly_parser.add_argument_group('assembly arguments')

assembly_group.add_argument('-min', metavar='INT', type=int, default=80, dest='min_overlap',
							help="minimun length of overlap [80]")

assembly_group.add_argument('-max', metavar='INT', type=int, default=90, dest='max_overlap',
							help="maximum length of overlap [90]")

assembly_group.add_argument('-oid', metavar='FLOAT', type=float, default=0.95, dest='overlap_identity',
							help="minimun identity of overlap region [0.95]")

assembly_group.add_argument('-tp', metavar='INT', type=int, dest='cluster_number_needKeep',
							help="how many clusters using in assembly. default=2")

assembly_group.add_argument('-ab', metavar='INT', type=int, dest='abundance_threshod',
							help="keep all clusters to assembly if its abundance >=INT ")

assembly_group.add_argument('-seqs_lim', metavar='INT', type=int,
					 		default=0, help="reads number limitation. [0]")

assembly_group.add_argument('-len', metavar='INT', type=int, default=400, dest='standard_length',
							help="standard reads length [400]")

assembly_group.add_argument('-mode', metavar='INT', type=int, choices=[1,2], default=1,
							help="modle 1 is to cluster and keep most [-tp] abundance clusters,\n" +\
							"or clusters abundance more than [-ab], and then make a consensus\n" +\
							"sequence for each cluster. modle 2 is directly to make only \n" +\
							"consensus sequence without clustering.")

assembly_group.add_argument('-rc', dest='reads_check', action="store_true",
							help="whether to check amino acid translation for reads")

assembly_group.add_argument('-cc', dest='coi_check', action="store_true",
							help="whether to check final COI contig's amino acid translation") 

assembly_group.add_argument('-codon', metavar='INT', type=int, dest="codon_table", default=5,
							help="codon table using to check translation [5], \n" +\
							"by the way, table [4,5] have same effect for COI gene.")

assembly_group.add_argument('-frame', metavar='INT', type=int,choices=[0,1,2], 
							default=1, help="translation start shift [1]") 


## only assembly need
only_assembly_parser = argparse.ArgumentParser(add_help=False)
only_assembly_group = only_assembly_parser.add_argument_group('when only run assembly arguments')

only_assembly_group.add_argument('-list', metavar='FILE', type=str, required=True,
							help="input file, fastq file list. [required]")
#------------------------------------------------------------------------------------------------

###############################################################################
#####----------------------- main subcommand parsers --------------------######

description = """

Description

	An automatic pipeline for HIFI-SE400 project, including filtering raw reads, 
	assigning reads to samples, assembly HIFI barcodes (COI sequences). 

Version
	1.0 2018-11-3

Author
	yangchentao at genomics.cn, BGI.
	mengguanliang at genomics.cn, BGI.
	
"""

parser = argparse.ArgumentParser(prog="HIFI-SE", description=description,
#								formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#								formatter_class=argparse.RawDescriptionHelpFormatter)
								formatter_class=argparse.RawTextHelpFormatter)

subparsers = parser.add_subparsers(dest='command')

########## subcommmands ###########

## all subcommand
parser_all = subparsers.add_parser("all", parents=[common_parser, index_parser, soft_parser, 
									filter_parser, assign_parser, assembly_parser], 
									formatter_class=argparse.RawTextHelpFormatter,
									help="run filter, assign and assembly")

## filter subcommand
parser_filter = subparsers.add_parser("filter", parents=[common_parser, filter_parser], 
									formatter_class=argparse.RawTextHelpFormatter,
									help="filter raw reads")

## assign subcommand
parser_assign = subparsers.add_parser("assign", parents=[common_parser,index_parser,
										only_assign_parser,assign_parser],
										formatter_class=argparse.RawTextHelpFormatter,
										help="assign reads to samples")

## assembly subcommand
parser_assembly = subparsers.add_parser("assembly", parents=[common_parser,index_parser,
										only_assembly_parser,soft_parser, assembly_parser],
										formatter_class=argparse.RawTextHelpFormatter, 
										help="do assembly from input fastq reads,\n" +\
										"output HIFI barcodes.")

## BOLD_identification
parser_bold = subparsers.add_parser("bold_identification", parents=[ ],
									formatter_class=argparse.RawTextHelpFormatter, 
										help="do taxa identification on BOLD system,\n")

###############################################################################
###############################################################################
#####---------------------- program execution start ----------------------#####

if len(sys.argv) == 1:
	parser.print_help()
	parser.exit()

#------------------------BOLD identification--------------------------#
if sys.argv[1] == 'bold_identification':
#if args.command == 'bold_identification':
    sys.argv = sys.argv[1:]
    sys.exit(bold_identification())

#----------------------- BOLD identification end ---------------------#

if sys.argv[1] in ['all', 'filter', 'assign', 'assembly']:
	args = parser.parse_args()
else:
	parser.print_help()
	parser.exit()

#-----------------------arguments checking-----------------------#
## softwares and databases
errors_found = 0
def check_program_involed(cmd):	
    result = subprocess.call('type %s' % cmd, shell = True, 
    		stdout = subprocess.PIPE, stderr = subprocess.PIPE) == 0
    if result:
    	return 0
    else:
    	print(cmd, " not found!", file=sys.stderr)
    	return 1

def files_exist_0_or_1(filelist):
	num = 0
	for file in filelist:
		if os.path.exists(file):
			num += 1
		else:
			print("%s doesn't exist!" % file, file=sys.stderr)
	if len(filelist) == num:
		return 0
	else:
		return 1

#----------------------------------------------------------------
## file existing check
if args.command == 'all':
	errors_found += files_exist_0_or_1([args.raw, args.primer])
elif args.command == 'filter':
	errors_found += files_exist_0_or_1([args.raw])
elif args.command == 'assign':
	errors_found += files_exist_0_or_1([args.primer])
elif args.command == "assembly":
	errors_found += files_exist_0_or_1([args.list])
else:
	parser.print_help()
	parser.exit()

if args.command in ['all', 'assembly']:
	vsearch = "vsearch"
	if hasattr(args, "vsearch"):
		if args.vsearch:
			vsearch = os.path.join(args.vsearch, vsearch)
	errors_found += check_program_involed(vsearch)

if errors_found > 0:
	parser.exit("Errors found! Exit!")


if args.outpre[-1:] == "/":
	print("outpre is in bad format!")
	eixt()

#-----------------------functions for filtering------------------#
def exp_e(q):

	exp = 0
	tmp = list(q)
	ascill = [ ord(n) - 33 for n in tmp ]

	for i in ascill:
		exp += 10 **(-i/10)

	return exp

def lowquality_rate(qual_str,cut_off):
	low_base = 0
	tmp = list(qual_str)
	ascill = [ ord(n) - 33 for n in tmp ]

	for i in ascill:
		if i < cut_off:
			low_base += 1

	low_rate = low_base / len(qual_str)
	return low_rate

#----------------------functions for assigning-------------------#

def comp_rev(sequence):
    # make a sequence complement #
    sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return sequence.upper()[::-1]

#----------------------functions for assembling------------------#

def complementation(sequence):
	# make a sequence complement #
	sequence.upper()
	sequence = sequence.replace('A', 't')
	sequence = sequence.replace('T', 'a')
	sequence = sequence.replace('C', 'g')
	sequence = sequence.replace('G', 'c')
	return sequence.upper()

def complement_and_reverse(reads_list):
	# make a sequence reverse and complement #
	
	new_reads_list = []
	for read in reads_list:
		read = complementation(read)
		new_reads_list.append(read[::-1])

	return new_reads_list

# ----count matched bases of two sequences----#
def match(str1,str2):
	matched = 0
	for base in range(len(str1)) :
		if str1[base] == str2[base]:
			matched += 1
	identity = matched / len(str1)
	return identity

#---------TranslateDNASeq------------#

def TranslateDNASeq(seq):
	l_dna = len(seq)
	if l_dna % 3 is not 0:
		seq = seq[:-(l_dna % 3)]
		#print("your sequence lenght is not tripple" + \
		#	"but no worries, I have trimmed well format")
	coding_dna = Seq(seq, generic_dna)
	protein = coding_dna.translate(table=args.codon_table)
	if "*" in protein:
		return False
	else:
		return True

#----------read fastq file-----------#
def read_fastq(fastq_file,ori):
	if os.path.splitext(fastq_file)[-1][1:] == "gz":
		fh_file =  gzip.open(fastq_file,'rt')
	else:
		fh_file = open(fastq_file,'r')

	reads_count = 0
	good = 0
	bad_length = 0
	seq_checked = []
	#fastq_err = "It's not a correct fastq format.\n"
	head = fh_file.readline().strip()
    # bug! id -> head
	while head:
		reads_count += 1
		sequence = fh_file.readline().strip()
		fh_file.readline()
		qual = fh_file.readline().strip()
		head = fh_file.readline().strip()
		seq_len = len(sequence)
		if  seq_len != args.standard_length :
			bad_length += 1
			continue

		#check translation
		if args.reads_check :
			tmp = sequence
			if ori == 'f':
				#forward primer length is 25
				tmp_remove = args.index + 25 + args.frame
				tmp = tmp[tmp_remove:]
				if TranslateDNASeq(tmp):
					seq_checked.append(sequence)
					good += 1
			else:
				#reverse primer length is 26
				#triming; complementation; triming end; reverse
			
				tmp_remove = args.index + 26 ;
				tmp = tmp[tmp_remove-1:]
				# complementation #
				tmp = complementation(tmp);
				# triming end #
				if seq_len % 3  is not 0 :
					tmp = tmp[0:-(seq_len % 3)]
				# reverse #
				tmp = tmp[::-1];
				
				if TranslateDNASeq(tmp):
					seq_checked.append(sequence)
					good += 1
		else:
			# do not check translation #
			seq_checked.append(sequence)
		
		if args.seqs_lim and  reads_count >= args.seqs_lim:
			break
	
	fh_file.close()
	if args.reads_check:
		fh_log.write(ori + ": total reads input:" + str(reads_count) +\
                    "\tcheck codon good:" + str(good) + "\t" +\
                    "with wrong length reads:" + str(bad_length) + "\n")
	else:
		fh_log.write(ori + ": total reads input:" + str(reads_count) +"\n")

	return seq_checked
#---------------COI_check------------#
def COI_check(contig):
	for_trim = args.index + 25 + 1;
	rev_trim = args.index + 26;
	contig = contig[for_trim:]
    #contig = contig[0:for_trim] #bug!
	contig = contig[:-rev_trim]
	if TranslateDNASeq(contig):
		return True
	else:
		return False

#---------------depth_table----------#
def depth_table(seqs):
	llen = len(seqs[0])
	m = 0
	consensus = ''
	while m < llen:
		depth ={}
		for align in seqs:
			if align[m] in depth.keys():
				depth[align[m]] += 1
			else:
				depth[align[m]] = 1

		if 'N' in depth.keys():
			del depth['N']

		#sorted_base = sorted(depth.iteritems(), key=lambda x: x[1], reverse=True)#
		#sorted_base = sorted(depth, key=depth.__getitem__,reverse=True)
		sorted_base = sorted(depth, key=lambda k: (depth[k],k),reverse=True)

		# for some case, different bases have same deepth #
	
		consensus += sorted_base[0]
		
		depth ={}
		m += 1
	return consensus

def report_depth(table,title,seq,read_len,step,ori):
	depth_sum = {}
	reports = []
	four_bases = ('A','T', 'C', 'G')
	reports = table[seq]
	if ori == 'f':
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
		for x in range(step,read_len):
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
	#in mode 1, using all codon checked reads to cluster with 100% identity, this is easy to achieve.#
	abu = {}
	for s in seqs:
		if s in abu.keys():
			abu[s] += 1
		else:
			abu[s] = 1
	seqs = []

	#sort and cluster with 100% identity#
	#sorted_seq = sorted(abu, key=abu.__getitem__,reverse=True)
	sorted_seq = sorted(abu, key=lambda k: (abu[k],k),reverse=True)

	most_abuns = {}

	#keep top 1#
	if abu[sorted_seq[0]] > 1:
		#keep top 2, if top1 has same abundance with top2's#
		if abu[sorted_seq[0]] == abu[sorted_seq[1]]:
			most_abuns[sorted_seq[0]] = abu[sorted_seq[0]]
			most_abuns[sorted_seq[1]] = abu[sorted_seq[1]]
			
			print('Alert: top1 and top2 have same number of sequences!\n')
			fh_log.write('Alert: top1 and top2 have same number of sequences!\n') 

		else:
			most_abuns[sorted_seq[0]] =abu[sorted_seq[0]]
	else:
		# all reads are unique, so that's wrong! #
		most_abuns[sorted_seq[0]] = abu[sorted_seq[0]]
		print("all reads in file are uniqe! clustering is meaningless!\n")

	abu = {}

	return most_abuns

def merge_matrix(matrix1,matrix2):

	for i in range(len(matrix2)):
		matrix1.append(matrix2[i])

	return matrix1

def mode_vsearch(seqs):
	# in mode 1, but clustering identity is not 100%, it is hard to archieve by perl, so I use VSEARCH to make it.
	# I just keep three sequences with the highest abundance after sorting all clusters.
	#
	# HASH table %cluster_tables contains pairs, which key = one cluster's consensus sequence, val = all sequences in
	# this cluster with 2D array format.
	cluster_tables = {} # 'consensus sequence' -> 'bases matrix build from same cluster'
	pid = os.getpid()
	temp_fasta = "temp.fa" + "." + str(pid)
	temp_uc  = "temp.uc" + "." + str(pid)
	with open(temp_fasta, 'w') as TM:
		for i in range(len(seqs)):
			TM.write(">" + str(i) + "\n" + seqs[i] + "\n")


	vsearch_cmd = vsearch + \
				" --cluster_fast " +  temp_fasta + \
				" --threads " + str(args.threads) + " --quiet " +\
				" --uc " + temp_uc +\
				" --id " + str(args.cluster_identity)
	#print("now run: " + vsearch_cmd)
	subprocess.call(vsearch_cmd, shell=True)

	clusters = {}
	count = {}

	# open temp.uc and statistic each cluster's abundance.
	with open(temp_uc,'r') as uc:
		for line in uc.readlines():
			if line[0] is not "H":  # there are "H","S","C" in the head of line
				continue
			array = line.split()
			if array[9] in clusters.keys():
				clusters[array[9]].append(array[8])
				count[array[9]] += 1
			else:
				clusters[array[9]] = []
				clusters[array[9]].append(array[9])
				count[array[9]] = 1

	#sorting clusters by abundance.#
	#sorted_clusters = sorted(count, key=count.__getitem__,reverse=True)
	sorted_clusters = sorted(count, key=lambda k:(count[k],k),reverse=True)


	if args.cluster_number_needKeep:
		#if set "-tp", keep top N clusters to assembly
		keep = args.cluster_number_needKeep
		sorted_clusters = sorted_clusters[0:keep]
	elif args.abundance_threshod:
		#if set "-ab", keep all clusters of abundance > ab
		while(sorted_clusters):
			item = sorted_clusters.pop()
			if count[item] < args.abundance_threshod:
				pass
			else:
				sorted_clusters.append(item)
				break
	else:
		# if set nothing, I will set it to -tp 2; 
		# if second most abundant sequence less than 1/10 of first, remove it!# of course it is just for when tp==2
		sorted_clusters = sorted_clusters[0:2]
		if count[sorted_clusters[1]] < count[sorted_clusters[0]]/10:
			sorted_clusters.pop()


	for  k in sorted_clusters:
		k_seqs = []
		matrix = []
		for s in clusters[k]:
			s = int(s)
			k_seqs.append(seqs[s])
			#split reads into bases-arrays#
			bases = list(seqs[s])
			matrix.append(bases)

		con = depth_table(k_seqs)

		#!!!---BUG---report
		# if make a consensus sequence for each cluster, you may face this problem: different clusters have same
		# consensus sequence, so previous key will be masked.
		# On! fuck!
		# 
		# Take it easy! I have a solution! why not to merge two pairs which have same keys.
		#
		if con in cluster_tables.keys():
			#merge this matric and previous same matrix#
			flash = merge_matrix(cluster_tables[con],matrix)
			cluster_tables[con] = flash
		else:
			cluster_tables[con] = matrix

	return cluster_tables

def mode_consensus(seqs):
	#In mode 2, this is only table.#

	matrix = []
	cons_table = {}
	for i in seqs:
		a = list(i)
		matrix.append(a)

	consensus = depth_table(seqs)
	cons_table[consensus] = matrix
	return cons_table;

#------------------------filter process--------------------------#
if args.command in ['all','filter']:

	filtered_outfile = args.outpre + "_filter_highQual.fq"
	if os.path.exists(filtered_outfile):
		print("WARRNING: " + filtered_outfile + " exists! now overwriting")
	out = open(filtered_outfile,'w')

	#Read sequences.
	err = open(args.outpre +  "_filter_lowQual.fastq",'w')
	log = open(args.outpre + "_filter_log.txt",'w')

	if os.path.splitext(args.raw)[-1][1:] == "gz":
		fh =  gzip.open(args.raw,'rt')
	else:
		fh = open(args.raw,'r')
	
	if args.expected_err and args.quality:
		print( "Bad arguments:\n\t" +\
				"-e argument is confilicting with -q," +\
				" can not using in the same time" )
		exit()
	elif args.quality:
		high_qual = args.quality[0]
		low_qual_cont = args.quality[1]/100
		filter_type = 2
		log.write("Filtering by quality score:\targs.quality\n")
	else:
		filter_type = 1
		log.write("Filtering by expected_err:{}".format(args.expected_err) + "\n")

	print("Filtering start: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	
	total = 0
	clean = 0
	nn = 0
	id = fh.readline().strip()
	if id[0] != '@':
		print("ERROR: {} is not a correct fastq format".format(args.raw))
		exit()
	while id:
		total += 1
		seq = fh.readline().strip()
		fh.readline().strip()
		qual = fh.readline().strip()
		N_count=seq.count('N')

		if N_count < args.n:
			if filter_type == 1:

				if exp_e(qual) <= args.expected_err :
					out.write(id + "\n" + seq + "\n" + "+\n" + qual + "\n" )
					clean += 1
				else:
					err.write(id + "\n" + seq + "\n" + "+\n" + qual + "\n")
			else:
				if lowquality_rate(qual,high_qual) > low_qual_cont :
					out.write(id + "\n" + seq + "\n" + "+\n" + qual + "\n" )
					clean += 1
				else:
					err.write(id + "\n" + seq + "\n" + "+\n" + qual + "\n")
		else:
			nn += 1

		id = fh.readline().strip()

	log.write("total reads:\t{}".format(total) + "\n")
	log.write("containing N reads:\t{}".format(nn) + "\n")
	log.write("clean reads:\t{}".format(clean))
	
	fh.close()
	err.close()
	log.close()
	out.close()
	
	print("Filtering done: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

#------------------------assign process--------------------------#

if args.command in ['all','assign']:
	print("Assigning start: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	if args.command == "all":
		args.fq = filtered_outfile
		assigned_outdir = os.path.abspath(args.outpre + "_assign")
		ErrFile = open(assigned_outdir + "_assign_err.fasta",'w')

	elif args.command == "assign":
		assigned_outdir = os.path.abspath(args.outdir)
		ErrFile = open(assigned_outdir + "_err.fasta",'w')

	indexlen = args.index

	if os.path.exists(assigned_outdir) == False :
	    os.mkdir(assigned_outdir)

	# check primer list number
	less_cmd = "less -S " + args.primer + "|wc -l"
	priwc = subprocess.check_output(less_cmd,shell=True)
	primer_lines = priwc.decode('utf-8')
	if int(primer_lines) % 2 != 0:
	    print("the primer file need to have each forward and reverse primer")
	    exit()

	pris = {}
	indp = {}
	FH = {}

	with open(args.primer,'r') as p:
	    for i in p.readlines():
	        i = i.strip()
	        arr = i.split()
	        if len(arr) != 2:
	            print("primer set is not well-formated")
	            exit()
	        sam = arr[0]
	        ipr = arr[1]
	        FH[ipr]= arr[0]
	        if sam in indp.keys():
	            print( arr[0] + "show twice in primer set")
	        else:
	            ori = sam[0:3]
	            num = sam[-3:]

	            pris[num] = {}
	            pris[num][ori] = ipr
	            indp[sam] = ipr
	            
	        if 'For' in sam:
	            plenf = len(ipr)
	            primerF = ipr[indexlen:]

	        if 'Rev' in sam:
	            plenr = len(ipr)
	            primerR = ipr[indexlen:]

	neg_priF = comp_rev(primerF)
	neg_priR = comp_rev(primerR)

	assigned_list = args.outpre + "_assignlist"

	with open(assigned_list,'w') as ls:
	    sorted_sample = sorted(pris.keys()) 
	    for s in sorted_sample:
	        ls.write(assigned_outdir + "/For" + s + ".fastq" + "\t" +\
	            	assigned_outdir + "/Rev" + s + ".fastq" + "\n")

	filehandle = {}
	for sam in indp.keys():
	    filehandle[sam] = open(assigned_outdir + "/" + sam + ".fastq",'w')

	seqnum=0
	err = 0
	assigned = 0
	count_assigned = {}
	count_total = {}

	with open(args.fq,'r') as fh:
	    head = fh.readline().strip()
	    while head:
	        seq = fh.readline().strip()
	        fh.readline().strip()
	        qual = fh.readline().strip()
	        head = fh.readline().strip()
	        seqnum += 1
	        headf = seq[0:plenf]
	        headr = seq[0:plenr]
	        if 'N' in seq:
	            continue
	        if headf in  FH:
	            tmp = seq[plenf:]
	            if FH[headf] in count_assigned:
	                count_total[FH[headf]] += 1
	            else:
	                count_total[FH[headf]] = 1
	            
	            if primerF not in tmp and primerR not in tmp and neg_priF not in tmp and neg_priR not in tmp:
	                
	                filehandle[FH[headf]].write("@" + FH[headf] + "_" + str(seqnum) + "\n" + seq + "\n")
	                filehandle[FH[headf]].write("+\n"  + qual + "\n")
	                
	                if FH[headf] in count_assigned:
	                    count_assigned[FH[headf]] += 1
	                else:
	                    count_assigned[FH[headf]] = 1
	                assigned += 1
	            else:
	                err += 1
	                ErrFile.write(">" + str(seqnum) + "\n" + seq + "\n")

	        elif headr in FH :
	            tmp = seq[plenr:]
	            if FH[headr] in count_assigned:
	                count_total[FH[headr]] += 1
	            else:
	                count_total[FH[headr]] = 1
	            
	            if primerF not in tmp and primerR not in tmp and neg_priF not in tmp and neg_priR not in tmp:
	                filehandle[FH[headr]].write("@" + FH[headr] + "_" + str(seqnum) + "\n" + seq + "\n")
	                filehandle[FH[headr]].write("+\n"  + qual + "\n")
	                if FH[headr] in count_assigned:
	                    count_assigned[FH[headr]] += 1
	                else:
	                    count_assigned[FH[headr]] = 1
	                assigned += 1
	            else:
	                err += 1
	                ErrFile.write(">" + str(seqnum) + "\n" + seq + "\n")
	            
	        else:
	            ErrFile.write(">" + str(seqnum) + "\n" + seq + "\n")
	            err += 1

	ErrFile.close()
	with open(args.outpre + ".assign.log",'w') as log:

	    log.write("total reads:\t{}\n".format(seqnum)) 
	    log.write("err reads:\t{}\n".format(err))
	    log.write("assigned:\t{}\n".format(assigned))

	    for i in sorted(count_assigned.keys()):
	        log.write(i + "\t" + str(count_total[i]) + "\t" + str(count_assigned[i]) + "\n")

	print("Assigning done: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

#------------------------assembly process--------------------------#
if args.command in ['all','assembly']:
	print("Assembling start: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	if args.cluster_number_needKeep and args.abundance_threshod:
		print( "Bad arguments:\n\t" +\
				"-tp argument is confilicting with -ab," +\
				" can not using in the same time" )
		exit()

	if args.min_overlap >83:
		print("For COI barcodes, by and large, overlaped length is 83 bp, so {} " +\
				"is not proper!".format(args.min_overlap))
		exit()
	if args.max_overlap < args.min_overlap:
		print("maximum overlap length must be large than minimun")
		exit()

	if args.command == 'all':
		args.list = assigned_list #list generated from assign step

	assembly_result = args.outpre + "_assembly.fa"
	if os.path.exists(assembly_result):
		print( assembly_result + " file exists!  overwriting ...")

	fh_out = open(assembly_result,'w')
	fh_log = open(args.outpre + "_assembly.log",'w')
	fh_depth = open((args.outpre + "_assembly.depth"),'w')

	if args.coi_check:
		out_checked = assembly_result + ".checked"
		fh_out_checked = open(out_checked,'w')

	if args.seqs_lim :
		fh_log.write("## reads input limitation: " + str(args.seqs_lim) + "\n")
	else:
		fh_log.write("## Using all reads to make consensus, no limitation\n")

	if args.reads_check:
		fh_log.write("## check codon translation = yes\n")
	else:
		fh_log.write("## check codon translation = no\n")

	fh_log.write("## consensus mode = " + str(args.mode) + "\n")

	if args.mode == 1:
	    fh_log.write("## clustering identity = " + str(args.cluster_identity) + "\n") 

	fh_log.write("## overlaping identity = " + str(args.overlap_identity) + "\n")
	fh_log.write("## min overlap = " + str(args.min_overlap) + "\n")
	fh_log.write("## max overlap = " + str(args.max_overlap) + "\n")

	list_format_info = "your list is not well formated!" +\
	                    "for example:\n\t" +\
						"/path/test_For001.fastq\t/path/test_Rev001.fastq"


	#--------------main-----------------------#
	barcodes_count = 0
	try:
		with open(args.list) as fh_list:
			lines = fh_list.readlines()
	except FileNotFoundError:
		print("can not find " + args.list)
		exit(0)


	for line in lines:
		line = line.rstrip()
		tmp_list = line.split()
		if len(tmp_list) is not 2:
			print(list_format_info)
			exit(0)
		else:
			forward = tmp_list[0]
			reverse = tmp_list[1]
			for_name = os.path.basename(forward).split('.')[0]
			rev_name = os.path.basename(reverse).split('.')[0]
			outname = for_name + "_" + rev_name
			fh_log.write("//processing " + outname + "\n")

			seq_checked_for = read_fastq(forward,'f')
			seq_checked_rev = read_fastq(reverse,'r')

			if len(seq_checked_for) == 0 or len(seq_checked_rev) ==0:
				fh_log.write("eithor Forward or Reverse file is empty!" + "\n")
				continue

			#here table_f and table_r are two quotes of two dicts, key is consensus sequence
			#value is a 2D array, including each cluster's bases.

			if args.mode == 1 and args.cluster_identity == 1:
				table_f = mode_identical(seq_checked_for)
				seq_checked_rev = complement_and_reverse(seq_checked_rev)
				table_r = mode_identical(seq_checked_rev)
			
			elif args.mode == 1 and args.cluster_identity < 1:
				table_f = mode_vsearch(seq_checked_for)
				seq_checked_rev = complement_and_reverse(seq_checked_rev)
				table_r = mode_vsearch(seq_checked_rev)

			else:
				#mod == 2
				#in mod2, though there is only one sequence, also using a array to store.
				table_f = mode_consensus(seq_checked_for)
				seq_checked_rev = complement_and_reverse(seq_checked_rev)
				table_r = mode_consensus(seq_checked_rev)

			#sort array by it's abundance#
			consensus_for = sorted(table_f.keys(), key=lambda x : len(table_f[x]), reverse=True)
			consensus_rev = sorted(table_r.keys(), key=lambda x : len(table_r[x]), reverse=True)

			##-----------anchoring overlap site--------#
		
			short_outname = outname[-3:]
			i = 0
			j = 0
			len_conF = len(consensus_for)
			len_conR = len(consensus_rev)
			for i in range(0,len_conF):
				cluster_f = consensus_for[i]
				pos1 = i +1
				abundance_f = len(table_f[cluster_f])
				
				fh_log.write(">" + short_outname + "_f_" + str(pos1) +"\tsize=" + str(abundance_f) + "\n" + cluster_f+"\n")

				for j in range(0,len_conR):
					cluster_r = consensus_rev[j]
					pos2 = j + 1
					abundance_r = len(table_r[cluster_r])

					# print just once #
					if pos1 == 1:
						fh_log.write(">" + short_outname + "_r_" + str(pos2) + "\tsize=" + str(abundance_r) + "\n" + cluster_r +"\n") 
					
					c_size = (abundance_f + abundance_r)/2
					c_size = str(c_size)
					
					read0 = cluster_f[-args.max_overlap:]
					read1 = cluster_r[0:args.max_overlap]
					
					singal = 0
					overlaps = {}
					for s in range(args.min_overlap, args.max_overlap + 1):
					
						l0 = read0[-s:]
						l1 = read1[0:s]
						tmp_identity = match(l0,l1)
						if tmp_identity == 1:
							overlaps[s] = 1
							fh_log.write( str(pos1) + "-" + str(pos2) + " ---> find position " +\
										str(s) +" with 100% mathch\n") 
							#find best result, so exit loop #
							break

						elif tmp_identity >= args.overlap_identity:
							
							overlaps[s] = tmp_identity 

					# find best overlaping result in all potenial positions 
					# candidates = sorted(overlaps.items(), lambda x, y: cmp(x[1], y[1]), reverse=True)#
					candidates = sorted(overlaps, key=overlaps.__getitem__,reverse=True)

					if len(candidates) > 0 :
						potenial = candidates[0]

						s0 = read0[-potenial:]
						s1 = read1[0:potenial]

						if  args.mode == 1 and args.cluster_identity== 1 :
							
							# if cid == 1, so depth of forward and reverse sequence is value of hash(table_f,table_r)
							#otherwise, mod==1,cid<1 or mod==2 have sample hash structure.

							fh_depth.write("# cid = 1, no need to report depth!") 

							if table_f[cluster_f] > table_r[cluster_r] :
								makeup_consensus = cluster_f + cluster_r[potenial-args.standard_length:]
							elif table_f[cluster_f] < table_r[cluster_r]:
								makeup_consensus = cluster_f[:args.standard_length-potenial] +  cluster_r
							else:
								makeup_consensus = cluster_f + cluster_r[potenial-args.standard_length:]
								fh_log.write("forward and reverse have same depth, use forward region\n") 
							
						else:
							correct = ''

							# report forward depth#
							title = short_outname + "_" + str(pos1) + "-" + str(pos2)
							report_depth(table_f, title, cluster_f,args.standard_length, potenial,'f')
							fh_depth.write(short_outname + "_" + str(pos1) + "-" + str(pos2) + "\t\t-----overlap start-----\n")

							##
							# compare each base from forward and reverse to keep one
							# forward == reverse
							# forward ne reverse, depth(forward) > depth(reverse)
							# forward ne reverse, depth(forward) < depth(reverse)
							
							for p in range(len(s0)):
							
								#site is changed, be careful!#
								tmp_loca0 = args.standard_length - potenial + p 

								# 2D array
								for_tab = table_f[cluster_f]
								rev_tab = table_r[cluster_r]
								(for_depth,rev_depth) = (0,0)
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
								fh_depth.write(short_outname + "_" + str(pos1) + "-" + str(pos2) + "\t" + str(tmp_loca1) +"\t") 
								four_bases = ('A','T','C','G')
								for base in four_bases:
									if base in sum.keys():
										fh_depth.write(base + ":" + str(sum[base]) + "\t")
									else:
										fh_depth.write(base+ ":0" + "\t") 
								# empty dict sum #
								sum = {}

								if s0[p] == s1[p] :
									correct += s0[p]
									info = s0[p] + "=" + s1[p]
								else:
									if for_depth > rev_depth:
										correct += s0[p]
										info = "[F=" + s0[p] +":"+ str(for_depth) + " > R=" + s1[p] +":" + str(rev_depth) + "]"
									elif for_depth == rev_depth:
										fh_log.write("short_outname\ttmp_loca1\tcould be a heterozygosis site!\n") 
										correct += s0[p]
										info = "[F=" + s0[p] +":"+ str(for_depth) + " = R=" + s1[p] +":" + str(rev_depth) + "]"
									else:
										correct += s1[p]
										info = "[F=" + s0[p] +":"+ str(for_depth) + " < R=" + s1[p] +":" + str(rev_depth) + "]"
							
								fh_depth.write(info + "\n") 

							fh_depth.write(short_outname + "_" + str(pos1) + "-" + str(pos2) + "\t\t-----overlap end-----\n") 
							# #report reverse depth #
							title = short_outname + "_" + str(pos1) + "-" + str(pos2) +"\t"
							report_depth(table_r, title, cluster_r, args.standard_length, potenial,'r')
							makeup_consensus = cluster_f[:args.standard_length-potenial] + correct + cluster_r[potenial-args.standard_length:]
						
						len_makeup_consensus = str(len(makeup_consensus))
						this_oid = overlaps[potenial] * 100
						this_oid = str("%.2f" % this_oid)

						# write into output file #
						barcodes_count += 1
						fh_out.write(">" + short_outname + "_" + str(pos1) + "-" + str(pos2) + ";" + str(abundance_f) + "_" +\
											str(abundance_r) + ";size=" + c_size + ";overPos=" + str(potenial) + ";oid=" + this_oid +\
											"%" + ";len=" + len_makeup_consensus + "\n" + makeup_consensus + "\n")
						# if check result, and ok so write into output_checked file #
						if args.coi_check and COI_check(makeup_consensus):
							fh_out_checked.write(">" + short_outname + "_" + str(pos1) + "-" + str(pos2) + ";" + str(abundance_f) +\
							"_" + str(abundance_r) + ";size=" + c_size + ";overPos=" + str(potenial) + ";oid=" +\
							this_oid + "%" + ";len=" + len_makeup_consensus + ";check=TRUE" + "\n" +\
							makeup_consensus + "\n")

					else:
						fh_log.write("<<" + str(pos1) + "-" + str(pos2) + " has no result!\n") 
	fh_out.close()
	if args.coi_check:
		fh_out_checked.close()

	if args.mode == 1 and args.cluster_identity < 1:
		rm_tmp_cmd = "rm temp.fa.* temp.uc.*"
		os.system(rm_tmp_cmd)

	print("Total barcodes generated: {}".format(barcodes_count))
	print("Assembling done: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

print("total run time: {}".format(time.time() -t))

