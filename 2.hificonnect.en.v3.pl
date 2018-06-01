#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::Seq;
use File::Basename;
use FindBin qw($Bin);
# 

=head1 Name

	hificonnect.pl

=head1 Description
	
	Due to connect HIFI barcode sequence by overlaping two consensus sequences which generated 
	from clustering method, in this program I use VSEARCH.  You can define the length of overlap,
	how many reads used to make clusters, and whether to check codon translation for PCG.

=head1 Version
	
	 Version: 0.3 ( May 23, 2018)
	 Author: yangchentao, yangchentao@genomics.cn

=head1 Usage 
	
	perl hifi.se400.pl [options] -list fq.list -mod 1

	--list   <str>   input file, fastq file list.
	--min    <int>   minimun length of overlap [60]
	--max    <int>   maximum length of verlap [90]
	--oid    <float> cutoff of identity of overlap region [0.85]
	--cid    <float> clustering identity rate [0.95]
	--tp     <int>   how many clusters using in assembly. [2]
	--limit  <y|n>   whether to limit sequencs number to save memory. [y|n]
	--seqs   <int>   reads number limitation. [50000]
	--index  <int>   index sequence lenght [4]
	--len    <int>   standard reads length [400]
	--mod    <1|2>   modle 1 is to cluster and keep 3 clusters with most abundance;
	                 modle 2 is to make a consensus sequence using all reads or all 
	                 codon checked reads if you set "-rc" ;
	                 * --cid is invaild in this mode;
	-out|o   <str>   output to file name [contig.fa]
	--rc      <?>    whether to check reads' codon translation [0]
	--cc      <?>    whether to check final COI contig's condon translation [1]
	--condon <int>   codon table using to check translation [5], by the way, table 4,5 have same effect for COI gene.
	--frame  <0|1|2> translation start shift [1]

=head1 Change logs

	1. remove mode 3;
	2. using vsearch program to cluster reads instead of hash table;
	3. using clusters depth info to recorrect consensus sequence instead of all reads;
	4. different length reads acceptable.

=cut

my ($help,$list,$min,$max,$limit,$reads_lim,$overlap_id,$cid,$mod,$len,$trim,$codon,$frame,$output,$toppick);
my ($read_check,$coi_check);

GetOptions(
	'help|h!' => \$help,
	'list=s' => \$list,
	'min:i'   => \$min,
	'max:i'   => \$max,
	'oid:f'    => \$overlap_id,
	'cid:f'   => \$cid,
	'limit|lim:s' => \$limit,
	'seqs:i'  => \$reads_lim,
	"len:i" => \$len,
	"index:i" => \$trim,
	"mod=i" => \$mod,
	"rc!" => \$read_check,
	"cc!" => \$coi_check,
	"tp:i" => \$toppick,
	"codon:i" => \$codon,
	"frame:i" => \$frame,
	"out|o:s"  => \$output
);


if (!$list or $help ){ 
	die `pod2text $0`;
}
# the region of overlap length
$min ||= 60;
$max ||= 90;
$mod ||= 1;
$codon ||= 5; # Invertebrate Mitochondrial
# COI 658bp barcode from second base to translate, now I didn't consider that condition indel occured before the start codon
# 1 + 219 *3 = 658
$frame ||= 1; 
$trim ||= 4; # barcode length. trim barocde before check codon.
$len ||= 400;
$overlap_id ||= 0.90;
$cid ||= 0.95;
$toppick ||= 2;
$output ||= "contig.fa";

if ($limit eq 'y'){
	$reads_lim ||= 50000;
	print "reads input limitation: $reads_lim\n";
}else{
	print "Using all reads to make consensus, no limitation\n";
}

if ($read_check) {
	print STDERR "check codon translation\n";
}else {
	print STDERR "do not check codon translation\n";
}

if ($coi_check) {
	open COI, ">$output.checked";
}
print STDERR "$output file exists!  overwriting ...\n" if (-e $output);

my $list_format_info = "your list is not well formated!
for example:

/path/test_For001.fastq\t/path/test_Rev001.fastq";

open OUT, ">$output";
open LOG, ">$output.log";
open LI,"$list" or die "No such file $list!";

print LOG "##consensus mode = $mod\n";
print LOG "##clustering identity = $cid\n";
print LOG "##overlaping identity = $overlap_id\n";

while (<LI>) {
	# read fastq file list, one sample per line
	chomp;
	my @a = split;
	my $for = $a[0];
	my $rev = $a[1];
	die "$list_format_info" unless ( @a == 2);
	my $c = $1 if (basename($for) =~ /(\S+)\./);
	my $d = $1 if (basename($rev) =~ /(\S+)\./);
	my $outname = "$c\_$d";
	print LOG "//processing $outname\n";

	my @seq_checked_for = &read_file($for,'f');

	my @seq_checked_rev = &read_file($rev,'r');

	# here $table_f and $table_r are two quotes of two HASH tables, key is consensus sequence
	# value is a 2D array, including each cluster's bases.
	
	my ($table_f,$table_r,@consensus_for,@consensus_rev);

	if ($mod == 1 && $cid == 1) 
	{
		$table_f = &mode_identical(@seq_checked_for);
		@seq_checked_rev = map {&comrev($_)} @seq_checked_rev;
		$table_r = &mode_identical(@seq_checked_rev);
		
	}elsif($mod == 1 && $cid < 1){

		$table_f = &mode_vsearch(@seq_checked_for);
		@seq_checked_rev = map {&comrev($_)} @seq_checked_rev;
		$table_r = &mode_vsearch(@seq_checked_rev);

	}else{
		# mod == 2

		# in mod2, though there is only one sequence, also using a array to store.
		
		$table_f = &mode_consensus(@seq_checked_for);
		@seq_checked_rev = map {&comrev($_)} @seq_checked_rev;
		$table_r = &mode_consensus(@seq_checked_rev);

	}

	# sort array by it's abundance
	@consensus_for = sort {$#{$$table_f{$b}} <=> $#{$$table_f{$a}}} keys %$table_f;
	@consensus_rev = sort {$#{$$table_r{$b}} <=> $#{$$table_r{$a}}} keys %$table_r;

	print STDERR "$outname\n";
	print STDERR "Making consensus done\n";
	print STDERR "Overlapping...\n";
	
	#-----------anchoring overlap site--------
	
	my $nammme = substr($outname,-3);
	for my $i(0..$#consensus_for){
		my $cluster_f = $consensus_for[$i];
		my $pos1 = $i +1;
		my $abundance_f = $#{$$table_f{$cluster_f}} + 1;
		print LOG ">$nammme\_f\_$pos1\tsize=$abundance_f\n$cluster_f\n";
		for my $j (0..$#consensus_rev){
			my $cluster_r = $consensus_rev[$j];
			my $pos2 = $j + 1;
			my $abundance_r = $#{$$table_r{$cluster_r}} + 1;
			print LOG">$nammme\_r\_$pos2\tsize=$abundance_r\n$cluster_r\n";
			my $c_size = ($abundance_f + $abundance_r) / 2;

			my $read0 = substr($cluster_f,-$max);
			my $read1 = substr($cluster_r,0,$max);

			my %smat;
			my %acceptable;
			my $singal = 0;
		TNT:for (my $s=$min; $s <= $max; $s++){
				my $l0 = substr($read0,-$s);
				my $l1 = substr($read1,0,$s);
				my $id = &match($l0,$l1);
				$smat{$s} = $id;
				if ($id == 1) {
					my $seq = &final_consensus($cluster_f,$cluster_r,$s);
					$singal = 1;
					my $len_seq = length $seq;
					if ($coi_check){
						if (COI_check($seq)){
							print COI ">$nammme\_$pos1-$pos2;$abundance_f\_$abundance_r;size=$c_size;overPos=$s;oid=100%;len=$len_seq;check=TRUE\n$seq\n";
							print OUT ">$nammme\_$pos1-$pos2;$abundance_f\_$abundance_r;size=$c_size;overPos=$s;oid=100%;len=$len_seq\n$seq\n";
						}
					}else{
						print OUT ">$nammme\_$pos1-$pos2;$abundance_f\_$abundance_r;size=$c_size;overPos=$s;oid=100%;len=$len_seq\n$seq\n";
					}

					print LOG "find position $s 100% mathch, jumping out loop!\n";
					last TNT; # find best result, so exit loop

				}elsif($id >= $overlap_id ){
					$acceptable{$s} = $id ;
				}
				
			}

			print STDERR "$pos1-$pos2 Finding overlap done\n";
			# find best overlaping result in all potenial positions
			my @candidates = sort {$acceptable{$b} <=> $acceptable{$a}} keys %acceptable;
			
			print LOG "<< $pos1-$pos2 has no result!\n" if ($singal == 0 && @candidates == 0 );

			if (@candidates > 0 && $singal == 0) {
				my $potenial = $candidates[0];
				my $correct;
				my $makeup_consensus;

				my $s0 = substr($read0,-$potenial);
				my $s1 = substr($read1,0,$potenial);

				if ($mod == 1 && $cid == 1){

					if ($$table_f{$cluster_f} > $$table_r{$cluster_r}) {
						$makeup_consensus = $cluster_f . substr($cluster_r,$potenial-$len);
					}elsif($$table_f{$cluster_f} < $$table_r{$cluster_r}){
						$makeup_consensus = substr($cluster_f,0,$len-$potenial) . $cluster_r;
					}else{
						$makeup_consensus = $cluster_f . substr($cluster_r,$potenial-$len);
						print LOG "forward and reverse have same depth, use forward region\n";
					}
				}else{
					my @res0 = split //, $s0;
					my @res1 = split //, $s1;
					my $count = 0;
					for my $p(0..$#res0){
						if ($res0[$p] eq $res1[$p]){
							$correct .= "$res0[$p]" ;
						}else {
							# site is changed, be careful!
							# !!!!
							my $tmp_loca0 = ($len - $potenial) + $p ;

#forward read
#								$q									
#0	1	2	3	4	5	6	7	8	9	10									forward length = 11
#A	T	C	G	G	A	A	G	T	T	A									reverse length =11
#							G	T	A	A	A	C	G	C	G	G	A		overlap = 4
#							0	1	2	3	4	5	6	7	8	9	10		$q = for_len - verlap + $p
#								$p											8 = 11 - 4 + 1
#reverse read


							# 2D array
							my @for_tab = @{$$table_f{$cluster_f}};
							my @rev_tab = @{$$table_r{$cluster_r}};
							my ($for_depth,$rev_depth) = (0,0);
							for my $x(0..$#for_tab){ 

								$for_depth ++ if ($for_tab[$x][$tmp_loca0] eq $res0[$p]);

							}

							for my $y(0..$#rev_tab){ 

								$rev_depth ++ if ($rev_tab[$y][$p] eq $res1[$p]); 

							}

							my $tmp_loca1 = $tmp_loca0 + 1;
							if ($for_depth > $rev_depth) {
								$correct .= $res0[$p];
							}elsif($for_depth == $rev_depth){
								print LOG "$nammme\t$tmp_loca1\tcould be a heterozygosis site!\n";
								$correct .= $res0[$p];
							}else{
								$correct .= $res1[$p];
							}
		
						}

					}
				}
				$makeup_consensus = substr($cluster_f,0,$len-$potenial) . $correct . substr($cluster_r,$potenial-$len);
				my $len_makeup_consensus = length $makeup_consensus;
				my $this_oid = $acceptable{$potenial} * 100;
				$this_oid = sprintf("%.2f",$this_oid);
				if ($coi_check){
					if (COI_check($makeup_consensus)){
						print OUT ">$nammme\_$pos1-$pos2;$abundance_f\_$abundance_r;size=$c_size;overPos=$potenial;oid=$this_oid%;len=$len_makeup_consensus\n$makeup_consensus\n";
						print COI ">$nammme\_$pos1-$pos2;$abundance_f\_$abundance_r;size=$c_size;overPos=$potenial;oid=$this_oid%;len=$len_makeup_consensus;check=TRUE\n$makeup_consensus\n";
					}
				}else{
					print OUT ">$nammme\_$pos1-$pos2;$abundance_f\_$abundance_r;size=$c_size;overPos=$potenial;oid=$this_oid%;len=$len_makeup_consensus\n$makeup_consensus\n";
				}
				
			}
		}
	}
}

close OUT;
close COI;
system("rm temp.fa.$$ temp.uc.$$");

#--------------------------------------------------------------------------

sub read_file{

	my $filename = shift;
	my $ori = shift;
	my @seq_checked;
	
		if ($filename =~ /.gz$/) 
	{
		open FQ,"gzip -dc $filename|" or die "Can not open $filename:$!\n";
	}else{
		open FQ,"$filename" or die "Can not open $filename:$!\n";
	}

	#Read sequences.
	my ($good,$count,$bad_length) = (0,0,0);
	while(my $t1=<FQ>) {
		$count ++;
		chomp($t1);
		if ($t1 !~ m/^@/)
		{
			print STDERR "ERROR:\n$filename: It's not a correct fastq format.\n";
			exit;
		}
		$t1 =~s/^@//;
		chomp(my $seq=<FQ>);
		<FQ>;
		chomp(my $qual=<FQ>);
		my $l = length $seq;
		my $ql = length $qual;
		
		if ($l != $len or $ql != $l) 
		{
			$bad_length ++;
			next;
		}

		my $tmp = $seq;

		# check translation
		if ($read_check){
		
			if ($ori eq 'f'){
				# forward primer length is 25
				my $tmp_remove = $trim + 25 + $frame;
				substr($tmp,0,$tmp_remove) = "";
				if (&TranslateDNASeq($tmp)) 
				{
					push @seq_checked, $seq;
					$good ++;
				}
				
			}else{
				# reverse primer length is 26
				# triming; complementation; triming end; reverse;
				my $tmp_remove = $trim + 26 ;
				# triming
				substr($tmp,0,$tmp_remove) = "";
				#complementation
				$tmp = &complementation($tmp);
				#triming end
				if ($l % 3  != 0){
					substr($tmp,-($l % 3)) = "";
				}
				# reverse
				$tmp = reverse ($tmp);
				
				if (&TranslateDNASeq($tmp))
				{
					push @seq_checked, $seq;
					$good ++;
				}				
			}
			
		}else{
			# do not check translation
			push @seq_checked, $seq;			
		}	
		
		last if ( $limit eq 'y' && $count >= $reads_lim);
	}
	close FQ;

	print LOG "$ori: total reads input:$count\tcheck codon good:$good\twith wrong length reads:$bad_length\n\n" if ($read_check);
	
	return @seq_checked;

}

#-----------


sub TranslateDNASeq
   {
	my $dna=shift;
	my $l_dna = length $dna;
	# to make sequence triple
	if ($l_dna % 3  != 0){
		substr($dna,-($l_dna % 3)) = "";
	}
	my $seqobj=Bio::Seq->new(-seq =>$dna, -alphabet =>'dna');
	my $prot_obj = $seqobj->translate( 
		                                 -codontable_id => $codon,
		                                 -terminator => 'U',
		                                 -unknown => '_',
		                                 -frame => 0
		                                 );
	my $pep = $prot_obj -> seq();

	($pep =~ /U/ || $pep =~ /_/) ? return 0 : return 1;
}

sub COI_check{
	my $contig = shift;
	my $for_trim = $trim + 25 + 1;
	my $rev_trim = $trim + 26;
	substr($contig,0,$for_trim) = "";
	substr($contig,-$rev_trim) = "";
	TranslateDNASeq($contig) ? return 1 : return 0;
}


sub depth_table
	{
		my @seqs = @_;
		my %depth;
		my $consensus;
		my $llen = length ($seqs[0]);
		for(my $y=0;$y<$llen;$y++){
			foreach my $align(@seqs){
				my $aa=substr($align,$y,1);
				$depth{$aa} ++; ## $a = (A|T|C|G|N);
			}

			delete $depth{"N"} if (exists $depth{'N'});

			my @sort_base = sort{$depth{$b} <=> $depth{$a} or $a cmp $b} keys %depth;

			# for some case, different bases have same deepth
		
			$consensus .= $sort_base[0];
			
			%depth = ();
			
		}
		
		return $consensus;
}

sub comrev {
	my $a = shift;
	chomp $a;
	$a =~ tr/ATCG/TAGC/;
	$a= reverse $a;
	return $a;
}

sub complementation {
	my $c = shift;
	chomp;
	$c =~ tr/ATCG/TAGC/;
	return $c;
}

sub match{
	my ($s0,$s1) = @_;
	my $len = length $s0;
	my @res0 = split //, $s0;
	my @res1 = split //, $s1;
	my $count = 0;
	for my $i(0..$#res0){
		$count ++ if ($res0[$i] eq $res1[$i]);
	}
	my $identity = $count / $len;
	return $identity;
}

sub final_consensus{
	my ($s0,$s1,$l) = @_;
	substr($s1,0,$l) = "";
	my $new = $s0.$s1;
	return $new;
}

sub mode_identical {
	
	# in mode 1, using all codon checked reads to cluster with 100% identity, this is easy to achieve.
	my @seqs = @_;
	my %abu;
	foreach (@seqs){
		$abu{$_} ++;
	}
	# sort and cluster with 100% identity
	my @sorted_seq = sort {$abu{$b} <=> $abu{$a} or $a cmp $b} keys %abu;
	my %most_abuns;

	if (@sorted_seq > 1) {
		if ($abu{$sorted_seq[0]} == $abu{$sorted_seq[1]}){
			
			$most_abuns{$sorted_seq[0]} = $abu{$sorted_seq[0]};
			$most_abuns{$sorted_seq[1]} = $abu{$sorted_seq[1]};
			
			print STDERR "top1 and top2 have same number of sequences!\n";
			print LOG "top1 and top2 have same number of sequences!\n";
			print LOG ">c1\tabundance=$abu{$sorted_seq[0]}\n$sorted_seq[0]\n";
			print LOG ">c2\tabundance=$abu{$sorted_seq[1]}\n$sorted_seq[1]\n";
		}else{
			$most_abuns{$sorted_seq[0]} = $abu{$sorted_seq[0]};
			print LOG ">c1\tabundance=$abu{$sorted_seq[0]}\n$sorted_seq[0]\n";
		}
		
	}else{
		$most_abuns{$sorted_seq[0]} = $abu{$sorted_seq[0]};
		print LOG ">c1\tabundance=$abu{$sorted_seq[0]}\n$sorted_seq[0]\n";
	}
	return \%most_abuns;
}

sub mode_vsearch {

	# in mode 1, but clustering identity is not 100%, it is hard to archieve by perl, so I use VSEARCH to make it.
	# I just keep three sequences with the highest abundance after sorting all clusters.
	#
	# HASH table %cluster_tables contains pairs, which key = one cluster's consensus sequence, val = all sequences in
	# this cluster with 2D array format.

	my @seqs = @_;
	my %cluster_tables;
	open TM,">temp.fa.$$" or die "$!";
	for my $i(0..$#seqs){
		print TM ">$i\n$seqs[$i]\n";
	}
	close TM;
	unless (-e "$Bin/vsearch"){
	
		die "can not find program: vsearch!\n";
	}
	system("$Bin/vsearch --cluster_smallmem temp.fa.$$ --threads 2 --quiet --uc temp.uc.$$ --id $cid");
	open UC, "< temp.uc.$$" or die "$!";
	my %clusters;
	my %count;
	while (<UC>) {
		chomp;
		next unless ($_=~/^H/);
		my @a = split /\t/;
		if (exists $clusters{$a[9]}) {
			push @{$clusters{$a[9]}}, $a[8];
			$count{$a[9]} ++;
		}else{
			push @{$clusters{$a[9]}}, $a[9];
			$count{$a[9]} = 1;
		}

	}
	close UC;

	# sorting clusters by abundance.
	my @sorted = sort {$count{$b} <=> $count{$a}} keys %count;
	
	# keep top two
	my $keep = $toppick - 1;
	@sorted = @sorted[0..$keep];

	# if second most abundant sequence less than 1/10 of first, remove it!
	if (@sorted == 2 && ($count{$sorted[1]} < ($count{$sorted[0]}/10)))
	{
		pop @sorted;
	}
	
	foreach my $k(@sorted) {
		
		my @k_seqs = ();
		my @matrix = ();
		foreach my $s(@{$clusters{$k}}){
			push @k_seqs, $seqs[$s];
			#split reads into bases-array
			my @bases = split //,$seqs[$s];
			push @matrix,[@bases];
		}

		my $con = &depth_table(@k_seqs);
	
		#!!!---BUG---report
		# if make a consensus sequence for each cluster, you may face this problem: different clusters have same
		# consensus sequence, so previous key will be masked.
		# On! fuck!
		# 
		# Take it easy! I have a solution! why not to merge two pairs which have same keys.
		if (exists $cluster_tables{$con}) {
			#merge
			my @flash = &merge_matrix(\@{$cluster_tables{$con}},\@matrix);
			@{$cluster_tables{$con}} = @flash;
		}else{
			@{$cluster_tables{$con}} = @matrix;	
		}
			
	}

	#traversal HASH table, write cluster information into LOG file.
	#my $id;
	#my @sort_by_abundance = sort {$#{$cluster_tables{$b}} <=>  $#{$cluster_tables{$a}} } keys %cluster_tables;
	#foreach my $key(@sort_by_abundance){
	#	$id ++;
	#	my @tt = @{$cluster_tables{$key}};
	#	my $this_cluster_abundance = @tt;
	#	print LOG ">cluster\_$id\tabundance=$this_cluster_abundance\n$key\n";

	#}
	return \%cluster_tables;
}

sub mode_consensus {
	# In mode 2, this is only table.
	my @seqs = @_;
	my $size = @seqs;
	my @matrix;
	my %cons_table;
	for my $i(0..$#seqs){
		my @a = split //,$seqs[$i];
		push @matrix,[@a];
	}
	my $consensus = &depth_table(@seqs);
	print LOG ">consensus\tsize=$size\n$consensus\n";
	@{$cons_table{$consensus}} = @matrix;
	return \%cons_table;
}

sub merge_matrix {
	my $matrix1 = shift;
	my $matrix2 = shift;
	for my $i(0..$#$matrix2){
		push @$matrix1, [@{$$matrix2[$i]}];
	}

	return @{$matrix1};
}
