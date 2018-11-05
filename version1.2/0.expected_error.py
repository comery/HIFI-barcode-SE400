import sys
import gzip
if len(sys.argv) < 3 :
	print("Usage: python " + sys.argv[0] + " [raw.fastq] " + "[expected error threshod]")
	exit()

filename = sys.argv[1]
threshod = int(sys.argv[2])
outfile = filename + ".clean.fq"


#Calculate phred score.
def exp_e(q):
	lenq = len(q)
	exp = 0
	phred = 0
	tmp = list(q)
	ascill = [ ord(n) - 33 for n in tmp ]

	for i in ascill:
		exp += 10 **(-i/10)

	return exp
#---------------------------------------------------------
#Read sequences.
err = open("err.fastq",'w')
log = open("log.txt",'w')
out = open(outfile,'w')

total = 0
nn = 0
clean = 0

with gzip.open(filename,'rt') as fh:
	id = fh.readline().strip()
	if id[0] != '@':
		print("ERROR: {} is not a correct fastq format".format(filename))
		exit()
	while id:
		seq = fh.readline().strip()
		fh.readline().strip()
		qual = fh.readline().strip()
		if 'N' in seq:
			err.write(id + "\n" + seq + "\n" + "+\n" + qual + "\n")
			nn += 1
	
		else:
			if exp_e(qual) < threshod :
				out.write(id + "\n" + seq + "\n" + "+\n" + qual + "\n" )
				clean += 1
		id = fh.readline().strip()

err.write("total reads:\t{}".format(total))
err.write("contain Ns reads:\t{}".format(nn))
err.write("clean reads:\t{}".format(clean))
err.close()
log.close()
out.close()


