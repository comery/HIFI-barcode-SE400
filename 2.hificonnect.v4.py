import os
import sys
import gzip
import argparse
import subprocess
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

'''parameters'''

parser = argparse.ArgumentParser(
	description="Due to connect HIFI barcode sequence by overlaping two " +\
	"consensus sequences which generated from clustering method, " +\
	"in this program I use VSEARCH.  You can define the length of overlap," +\
	"how many reads used to make clusters, and whether to check codon translation for PCG.")

parser.add_argument('--list', metavar='FILE', type=str,required=True,
					help="input file, fastq file list.")

parser.add_argument('--min', metavar='INT', type=int, default=60, dest='min_overlap',
					help="minimun length of overlap [60]")

parser.add_argument('--max', metavar='INT', type=int, default=90, dest='max_overlap',
					help="maximum length of overlap [90]")

parser.add_argument('--oid', metavar='FLOAT', type=float, default=0.85, dest='overlap_identity',
					help="cutoff of identity of overlap region [0.85]")

parser.add_argument('--cid', metavar='FLOAT', type=float, default=0.95, dest='cluster_identity',
					help="clustering identity rate [0.95]")

parser.add_argument('--tp', metavar='INT', type=int, default=2, dest='cluster_number_needKeep',
					help="how many clusters using in assembly. [2]")

parser.add_argument('--seqs_lim', metavar='INT', type=int,
					 default=50000, help="reads number limitation. [50000]")

parser.add_argument('--index', metavar='INT', type=int, default=4,
					help="index sequence lenght [4]")

parser.add_argument('--len', metavar='INT', type=int, default=400, dest='standard_length',
					help="standard reads length [400]")

parser.add_argument('--mode', metavar='INT', type=int, choices=[1,2], default=1,
					help="modle 1 is to cluster and keep 3 clusters with most abundance" +\
                         "modle 2 is to make a consensus sequence using all reads or all" +\
                         "codon checked reads if you set \"-rc\" " +\
                         "* --cid is invaild in this mode")

parser.add_argument('--out', metavar='FILE', type=str, default="contig.fa",
					help="output to file name [contig.fa]")


parser.add_argument('--rc', dest='reads_check', help="whether to check reads' codon translation",
					action="store_true")

parser.add_argument('--cc', dest='coi_check', help="whether to check final COI contig's condon translation",
					action="store_false") 

parser.add_argument('--codon', metavar='INT', type=int, dest="codon_table", default=5,
					help="codon table using to check translation [5], " +\
					"by the way, table 4,5 have same effect for COI gene.")

parser.add_argument('--frame', metavar='INT', type=int,choices=[0,1,2], default=1, help="translation start shift [1]") 

args = parser.parse_args()

'''------------------------main--------------------------'''

log_file = args.out+ ".log"
fh_log = open(log_file,'w')

if args.seqs_lim :
	print("## reads input limitation: " + str(args.seqs_lim), file=fh_log)
else:
	print("## Using all reads to make consensus, no limitation", file=fh_log)

if args.reads_check:
	print("## check codon translation = yes", file=fh_log)
else:
	print("## check codon translation = no", file=fh_log)

print("## consensus mode = " + str(args.mode), file = fh_log)
print("## clustering identity = " + str(args.cluster_identity), file = fh_log)
print("## overlaping identity = " + str(args.overlap_identity), file = fh_log)
print("## min overlap = " + str(args.min_overlap), file = fh_log)
print("## max overlap = " + str(args.max_overlap), file = fh_log)

if args.coi_check:
	out_checked = args.out + ".checked"
	fh_out_checked = open(out_checked,'w')
	if os.path.isfile(out_checked):
		print( str(args.out) + " file exists!  overwriting ...")

fh_out = open(args.out,'w')
fh_depth = open((args.out + ".depth"),'w')

list_format_info = "your list is not well formated!" +\
					"for example" +\
					"/path/test_For001.fastq\t/path/test_Rev001.fastq"

#if os.path.isfile(args.list):


def complementation(sequence):
	''' make a sequence complement '''
	sequence.upper()
	sequence = sequence.replace('A', 't')
	sequence = sequence.replace('T', 'a')
	sequence = sequence.replace('C', 'g')
	sequence = sequence.replace('G', 'c')
	return sequence.upper()


def complement_and_reverse(reads_list):
	''' make a sequence reverse and complement '''
	
	new_reads_list = []
	for read in reads_list:
		read = complementation(read)
		new_reads_list.append(read[::-1])

	return new_reads_list

''' ----count matched bases of two sequences----'''
def match(str1,str2):
	matched = 0

	for base in range(0,len(str1)) :
		if str1[base] == str2[base]:
			matched += 1
	identity = matched / len(str1)
	return identity


'''---------------TranslateDNASeq--------------------------'''

def TranslateDNASeq(seq):
	l_dna = len(seq)
	if l_dna % 3 != 0:
		seq = seq[:-(l_dna % 3)]
		'''print("your sequence lenght is not tripple" + \
			"but no worries, I have trimmed well format")'''
	coding_dna = Seq(seq, generic_dna)
	protein = coding_dna.translate(table=args.codon_table)
	if "*" in protein:
		return False
	else:
		return True


'''----------read fastq file-----------'''
def read_fastq(fastq_file,ori):
	if os.path.splitext(fastq_file)[-1][1:] == "gz":
		fh_file =  gzip.GzipFile(fastq_file)
	else:
		fh_file = open(fastq_file)

	reads_count = 0
	good = 0
	bad_length = 0
	seq_checked = []
	fastq_err = "It's not a correct fastq format.\n"
	head = fh_file.readline().strip()
	while id:
		reads_count += 1
		sequence = fh_file.readline().strip()
		fh_file.readline().strip()
		qual = fh_file.readline().strip()
		head = fh_file.readline().strip()
		if head == "":
			break
		
		if head[0] != "@" :
			print(fastq_err)
		
		seq_len = len(sequence)
		qual_len = len(qual)
		if  seq_len != 400 or qual_len != 400:
			bad_length += 1
			continue

		tmp = sequence

		''' check translation '''
		if args.reads_check :
			if ori == 'f':
				''' forward primer length is 25 '''
				tmp_remove = args.index + 25 + args.frame
				tmp = tmp[tmp_remove:]
				if TranslateDNASeq(tmp):
					seq_checked.append(sequence)
					good += 1
			else:
				'''
				reverse primer length is 26
				triming; complementation; triming end; reverse
				'''
				tmp_remove = args.index + 26 ;
				tmp = tmp[tmp_remove-1:]
				''' complementation '''
				tmp = complementation(tmp);
				''' triming end '''
				if seq_len % 3  != 0 :
					tmp = tmp[0:-(seq_len % 3)]
				''' reverse '''
				tmp = tmp[::-1];
				
				if TranslateDNASeq(tmp):
					seq_checked.append(sequence)
					good += 1
		else:
			''' do not check translation '''
			seq_checked.append(sequence)			
		
		if args.seqs_lim and  reads_count > args.seqs_lim:
			break
	
	fh_file.close()
	if args.reads_check:
		fh_log.write(ori + ": total reads input:" + str(reads_count) + "\tcheck codon good:" + str(good) + "\t" +\
					"with wrong length reads:" + str(bad_length) + "\n")
	
	return seq_checked


'''---------------COI_check--------------------------'''
def COI_check(contig):
	for_trim = args.index + 25 + 1;
	rev_trim = args.index + 26;
	contig = contig[for_trim:]
	contig = contig[0:for_trim]
	contig = contig[:-rev_trim]
	if TranslateDNASeq(contig):
		return True
	else:
		return False

'''---------------depth_table--------------------------'''
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

		'''sorted_base = sorted(depth.iteritems(), key=lambda x: x[1], reverse=True)'''
		sorted_base = sorted(depth, key=depth.__getitem__,reverse=True)

		''' for some case, different bases have same deepth '''
	
		consensus += sorted_base[0]
		
		depth ={}
		m += 1
	return consensus

def report_depth(table,title,seq,read_len, step,ori):
	depth_sum = {}
	reports = []
	four_bases = ('A','T', 'C', 'G')
	reports = table[seq]
	if ori == 'f':
		for_rest_len = read_len - step
		for x in range(0,for_rest_len):
			real_position = x + 1
			fh_depth.write(title + "\t" + str(real_position) + "\t") 
			for y in range(0,len(reports)):
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
		rev_rest_len = read_len -1
		for x in range(step,rev_rest_len):
			real_position = read_len - step + x + 1
			fh_depth.write(title + "\t" + str(real_position) + "\t")
			for y in range(0,len(reports)):
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
	
	'''in mode 1, using all codon checked reads to cluster with 100% identity, this is easy to achieve.'''
	abu = {}
	for s in seqs:
		if s in abu.keys():
			abu[s] += 1
		else:
			abu[s] = 1
	seqs = []

	'''sort and cluster with 100% identity'''
	sorted_seq = sorted(abu, key=abu.__getitem__,reverse=True)

	most_abuns = {}

	'''keep top 1'''
	print(abu[sorted_seq[0]])
	if abu[sorted_seq[0]] > 1:
		'''keep top 2, if top1 has same abundance with top2's'''
		if abu[sorted_seq[0]] == abu[sorted_seq[1]]:
			most_abuns[sorted_seq[0]] = abu[sorted_seq[0]]
			most_abuns[sorted_seq[1]] = abu[sorted_seq[1]]
			
			print('Alert: top1 and top2 have same number of sequences!\n')
			fh_log.write('Alert: top1 and top2 have same number of sequences!\n') 

		else:
			most_abuns[sorted_seq[0]] =abu[sorted_seq[0]]
	else:
		''' all reads are unique, so that's wrong! '''
		most_abuns[sorted_seq[0]] = abu[sorted_seq[0]]
		print("all reads in file are uniqe! clustering is meaningless!\n")

	abu = {}

	return most_abuns

def mode_vsearch(seqs):
	'''
	# in mode 1, but clustering identity is not 100%, it is hard to archieve by perl, so I use VSEARCH to make it.
	# I just keep three sequences with the highest abundance after sorting all clusters.
	#
	# HASH table %cluster_tables contains pairs, which key = one cluster's consensus sequence, val = all sequences in
	# this cluster with 2D array format.
	'''
	cluster_tables = {}
	pid = os.getpid()
	temp_fasta = "temp.fa" + "." + str(pid)
	temp_uc  = "temp.uc" + "." + str(pid)
	with open(temp_fasta, 'w') as TM:
		for i in range(0,len(seqs)):
			TM.write(">" + str(i) + "\n" + seqs[i] + "\n")

	vsearch_bin = os.path.join(sys.path[0],"vsearch")

	vsearch_cmd = vsearch_bin + \
				" --cluster_smallmem " +  temp_fasta + \
				" --threads 2 --quiet " +\
				" --uc " + temp_uc +\
				" --id 0.97 " 
	print("now run: " + vsearch_cmd)
	subprocess.call(vsearch_cmd, shell=True)

	clusters = {}
	count = {}

	with open(temp_uc,'r') as uc:
		for line in uc.readlines():
			if line[0] != "H":
				continue
			array = line.split()
			if array[9] in clusters.keys():
				clusters[array[9]].append(array[8])
				count[array[9]] += 1
			else:
				clusters[array[9]] = []
				clusters[array[9]].append(array[9])
				count[array[9]] = 1


	'''sorting clusters by abundance.'''
	sorted_clusters = sorted(count, key=count.__getitem__,reverse=True)
	
	'''keep top two'''
	keep = args.cluster_number_needKeep
	sorted_clusters = sorted_clusters[0:keep]

	'''if second most abundant sequence less than 1/10 of first, remove it!'''
	if len(sorted_clusters) == 2  and count[sorted_clusters[1]] < count[sorted_clusters[0]]/10:
		sorted_clusters.pop()

	for  k in sorted_clusters:
		k_seqs = []
		matrix = []
		for s in range(0,len(clusters[k])):
			k_seqs.append(seqs[s])
			'''split reads into bases-arrays'''
			bases = list(seqs[s])
			matrix.append(bases)

		con = depth_table(k_seqs)
		'''
		#!!!---BUG---report
		# if make a consensus sequence for each cluster, you may face this problem: different clusters have same
		# consensus sequence, so previous key will be masked.
		# On! fuck!
		# 
		# Take it easy! I have a solution! why not to merge two pairs which have same keys.
		'''
		if con in cluster_tables.keys():
			'''merge'''
			flash = merge_matrix(cluster_tables[con],matrix)
			cluster_tables[con] = flash
		else:
			cluster_tables[con] = matrix

	return cluster_tables

def mode_consensus(seqs):
	'''In mode 2, this is only table.'''

	matrix = []
	cons_table = {}
	for i in seqs:
		a = list(i)
		matrix.append(a)

	consensus = depth_table(seqs)
	cons_table[consensus] = matrix
	return cons_table;

def merge_matrix(matrix1,matrix2):
	
	for i in range(0,len(matrix2)):
		matrix1.append(matrix2[i])
	
	return matrix1


'''--------------main-----------------------'''
try:
	with open(args.list) as fh_list:
		lines = fh_list.readlines()
except FileNotFoundError:
	print("can not find " + args.list)

for line in lines:
	line = line.rstrip()
	tmp_list = line.split()
	if len(tmp_list) != 2:
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

		''' 
		here table_f and table_r are two quotes of two dicts, key is consensus sequence
		value is a 2D array, including each cluster's bases.
		'''

		if args.mode == 1 and args.cluster_identity == 1:
			table_f = mode_identical(seq_checked_for)
			seq_checked_rev = complement_and_reverse(seq_checked_rev)
			table_r = mode_identical(seq_checked_rev)
		
		elif args.mode == 1 and args.cluster_identity < 1:
			table_f = mode_vsearch(seq_checked_for)
			seq_checked_rev = complement_and_reverse(seq_checked_rev)
			table_r = mode_vsearch(seq_checked_rev)

		else:
		
			'''
			mod == 2
			in mod2, though there is only one sequence, also using a array to store.
			'''
			table_f = mode_consensus(seq_checked_for)
			seq_checked_rev = complement_and_reverse(seq_checked_rev)
			table_r = mode_consensus(seq_checked_rev)

		'''sort array by it's abundance'''
		consensus_for = sorted(table_f.keys(), key=lambda x : len(table_f[x]), reverse=True)
		consensus_rev = sorted(table_r.keys(), key=lambda x : len(table_r[x]), reverse=True)


		'''#-----------anchoring overlap site--------'''
	
		short_outname = outname[-3:]
		i = 0
		j = 0
		for i in range(0,len(consensus_for)):
			cluster_f = consensus_for[i]
			pos1 = i +1
			abundance_f = len(table_f[cluster_f])
			
			fh_log.write(">" + short_outname + "_f_" + str(pos1) +"\tsize=" + str(abundance_f) + "\n" + cluster_f+"\n")

			for j in range(0,len(consensus_rev)):
				cluster_r = consensus_rev[j]
				pos2 = j + 1
				abundance_r = len(table_r[cluster_r])

				''' print just once '''
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
						fh_log.write("---> find position " + str(s) +" with 100% mathch, jumping out loop!\n") 
						'''find best result, so exit loop '''
						break

					elif tmp_identity >= args.overlap_identity:
						
						overlaps[s] = tmp_identity 


				''' find best overlaping result in all potenial positions 
				candidates = sorted(overlaps.items(), lambda x, y: cmp(x[1], y[1]), reverse=True)'''

				candidates = sorted(overlaps, key=overlaps.__getitem__,reverse=True)

				if len(candidates) > 0 :
					potenial = candidates[0]

					s0 = read0[-potenial:]
					s1 = read1[0:potenial]

					if  args.mode == 1 and args.cluster_identity== 1 :
						'''
						if cid == 1, so depth of forward and reverse sequence is value of hash(table_f,table_r)
						otherwise, mod==1,cid<1 or mod==2 have sample hash structure.
						'''
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

						''' report forward depth'''
						title = short_outname + "_" + str(pos1) + "-" + str(pos2)
						report_depth(table_f, title, cluster_f,args.standard_length, potenial,'f')
						fh_depth.write(short_outname + "_" + str(pos1) + "-" + str(pos2) + "\t\t-----overlap start-----\n")

						'''
						compare each base from forward and reverse to keep one
						forward == reverse
						forward ne reverse, depth(forward) > depth(reverse)
						forward ne reverse, depth(forward) < depth(reverse)
						'''

						for p in range(0,len(s0)):
						
							'''site is changed, be careful!'''
							tmp_loca0 = args.standard_length - potenial + p 

							# 2D array
							for_tab = table_f[cluster_f]
							rev_tab = table_r[cluster_r]
							(for_depth,rev_depth) = (0,0)
							sum = {}

							sum[s0[p]] = 0
							sum[s1[p]] = 0
							for x in range(0,len(for_tab)):
								if for_tab[x][tmp_loca0] == s0[p]:
									for_depth += 1
									sum[s0[p]] += 1

							for y in range(0,len(rev_tab)):
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
							''' empty dict sum '''
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
						''' #report reverse depth '''
						title = short_outname + "_" + str(pos1) + "-" + str(pos2) +"\t"
						report_depth(table_r, title, cluster_r, args.standard_length, potenial,'r')
						makeup_consensus = cluster_f[:args.standard_length-potenial] + correct + cluster_r[potenial-args.standard_length:]
					
					len_makeup_consensus = str(len(makeup_consensus))
					this_oid = overlaps[potenial] * 100
					this_oid = str("%.2f" % this_oid)
					if args.coi_check:
						if COI_check(makeup_consensus):
							fh_out.write(">" + short_outname + "_" + str(pos1) + "-" + str(pos2) + ";" + str(abundance_f) + "_" +\
										str(abundance_r) + ";size=" + c_size + ";overPos=" + str(potenial) + ";oid=" + this_oid +\
										"%" + ";len=" + len_makeup_consensus + "\n" + makeup_consensus + "\n")
							fh_out_checked.write(">" + short_outname + "_" + str(pos1) + "-" + str(pos2) + ";" + str(abundance_f) +\
										"_" + str(abundance_r) + ";size=" + c_size + ";overPos=" + str(potenial) + ";oid=" +\
										this_oid + "%" + ";len=" + len_makeup_consensus + ";check=TRUE" + "\n" +\
										makeup_consensus + "\n")
					else:
						fh_out.write(">" + short_outname + "_" + str(pos1) + "-" + str(pos2) + ";" + str(abundance_f) + "_" +\
									str(abundance_r) + ";size=" + c_size + ";overPos=" + str(potenial) + ";oid=" + this_oid +\
									"%" + ";len=" + len_makeup_consensus + "\n" + makeup_consensus + "\n")

					
				else:
					fh_log.write("<<" + str(pos1) + "-" + str(pos2) + " has no result!\n") 
	fh_out.close()
	if args.coi_check:
		fh_out_checked.close()


	if args.mode == 1 and args.cluster_identity < 1:
		rm_tmp_cmd = "rm temp.fa.* temp.uc.*"
		os.system(rm_tmp_cmd)

