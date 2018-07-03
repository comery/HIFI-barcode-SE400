#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);


=head1 Name

	hificonnect.pl

=head1 Description
	
	DNA -> protein -> find most abundant seed  and remove outliers -> call consensus

=head1 Version
	
	 Version: 0.5 ( Jule 10, 2018)
	 Author: yangchentao, yangchentao@genomics.cn

=head1 Usage 
	
	perl hifi.se400.pl [options] -list fq.list 

	--list     <str>   input file, fastq file list.
	--min      <int>   minimun length of overlap [60]
	--max      <int>   maximum length of verlap [90]
	--oid      <float> cutoff of identity of overlap region [0.85]
	--cid      <float> clustering identity rate [0.95]
	--tp       <int>   how many clusters using in assembly. [2]
	--limit    <int>   limit sequencs number to save memory. 
	--index|i  <int>   index sequence lenght [4]
	--len      <int>   standard reads length [400]
	-out|o     <str>   output to file name [contig.fa]
	--rc       <?>    whether to check reads' codon translation [0]
	--cc       <?>    whether to check final COI contig's condon translation [1]
	--condon   <int>   codon table using to check translation [5], by the way, table 4,5 have same effect for COI gene.
	--frame    <0|1|2> translation start shift [1]

=head1 Change logs

	1. remove mode 3;
	2. using vsearch program to cluster reads instead of hash table;
	3. using clusters depth info to recorrect consensus sequence instead of all reads;
	4. different length reads acceptable.

=cut

my ($help,$list,$min,$max,$reads_lim,$overlap_id,$cid,$len,$trim,$codon,$frame,$output,$coi_check);

GetOptions(
	'help|h!' => \$help,
	'list=s' => \$list,
	'min:i'   => \$min,
	'max:i'   => \$max,
	'oid:f'    => \$overlap_id,
	'cid:f'   => \$cid,
	'limit|lim:i' => \$reads_lim,
	'len:i' => \$len,
	'index|i:i' => \$trim,
	'cc!' => \$coi_check,
	'codon:i' => \$codon,
	'frame:i' => \$frame,
	'out|o:s'  => \$output
);


if (!$list or $help ){ 
	die `pod2text $0`;
}
# the region of overlap length
$min ||= 60;
$max ||= 90;
$codon ||= 5; # Invertebrate Mitochondrial
# COI 658bp barcode from second base to translate, now I didn't consider that condition indel occured before the start codon
# 1 + 219 *3 = 658
$frame ||= 1; 
$trim ||= 4; # barcode length. trim barocde before check codon.
$len ||= 400;
$overlap_id ||= 0.90;
$output ||= "contig.fa";

if ($reads_lim){
	print "reads input limitation: $reads_lim\n";
}else{
	print "Using all reads to make consensus, no limitation\n";
}

if ($coi_check) {
	open COI, ">$output.checked";
}
print STDERR "$output file exists!  overwriting ...\n" if (-e $output);

my $list_format_info = "your list is not well formated!
for example:

/path/test_For001.fastq\t/path/test_Rev001.fastq";

system(" nohup sh $Bin/cpu_mem.record.sh $$ >> cpu_mem.stat$$ ");

open OUT, ">$output";
open LOG, ">$output.log";
open HET, ">$output.depth";
open LI,"$list" or die "No such file $list!";

print LOG "##overlaping identity = $overlap_id\n\n";

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

	my $prot_cluster_for = &read_file($for,'f');

	my $prot_cluster_rev = &read_file($rev,'r');


	my @normal_read_f = @{&remove_outliner($prot_cluster_for)};
	my @normal_read_r = @{&remove_outliner($prot_cluster_rev)};

	# here $table_f and $table_r are two quotes of two HASH tables, key is consensus sequence
	# value is a 2D array, including each cluster's bases.
	
	my ($table_f,$table_r,@consensus_for,@consensus_rev);


	$table_f = &mode_consensus(@normal_read_f);
	@normal_read_r = map {&comrev($_)} @normal_read_r;
	$table_r = &mode_consensus(@normal_read_r);

	@consensus_for = &depth_table(@normal_read_f);
	@consensus_rev = &depth_table(@normal_read_r);

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
			# print just once
			print LOG">$nammme\_r\_$pos2\tsize=$abundance_r\n$cluster_r\n" if ($pos1 == 1);
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

					#report forward depth
					&report_depth($table_f,"$outname\_$pos1-$pos2\t",$cluster_f,$s,'f');
					print HET "$outname\_$pos1-$pos2\t//\n";
					&report_depth($table_r,"$outname\_$pos1-$pos2\t",$cluster_r,$s,'r');

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

					print LOG "---> find position $s 100% mathch, jumping out loop!\n";
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

				#report forward depth
				&report_depth($table_f,"$outname\_$pos1-$pos2\t",$cluster_f,$potenial,'f');

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
						my %sum;

						for my $x(0..$#for_tab){ 

							$for_depth ++ if ($for_tab[$x][$tmp_loca0] eq $res0[$p]);
							$sum{$for_tab[$x][$tmp_loca0]} ++ ;

						}

						for my $y(0..$#rev_tab){ 

							$rev_depth ++ if ($rev_tab[$y][$p] eq $res1[$p]);
							$sum{$rev_tab[$y][$p]}++;

						}
						my $tmp_loca1 = $tmp_loca0 + 1;
						print HET "$outname\_$pos1-$pos2\t$tmp_loca1\t";
						my @four_bases = qw(A T C G);
						foreach my $base (@four_bases){
							
							exists $sum{$base} ? print HET "$base:$sum{$base}\t" : print HET "$base:0\t";
						}

						%sum = ();
						
						my $info;
						if ($for_depth > $rev_depth) {
							$correct .= $res0[$p];
							$info = "[F=$res0[$p]:$for_depth > R=$res1[$p]:$rev_depth]";
						}elsif($for_depth == $rev_depth){
							print LOG "$nammme\t$tmp_loca1\tcould be a heterozygosis site!\n";
							$correct .= $res0[$p];
							$info = "[F=$res0[$p]:$for_depth = R=$res1[$p]:$rev_depth]";
						}else{
							$correct .= $res1[$p];
							$info = "[F=$res0[$p]:$for_depth < R=$res1[$p]:$rev_depth]";
						}
						print HET "$info\n";
	
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

				#report reverse depth
				&report_depth($table_r,"$outname\_$pos1-$pos2\t",$cluster_r,$potenial,'r');
			}
			
		}
	}
}

close OUT;
close COI;

#--------------------------------------------------------------------------

sub read_file{

	my $filename = shift;
	my $ori = shift;
	my %prot_cluster;
	
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
		}else{
			$good ++;
		}

		my $pep = &TranslateDNASeq($seq);
			
		push @{$prot_cluster{$pep}}, $seq;
			

		last if ( $reads_lim && $count >= $reads_lim);
	}
	close FQ;

	print LOG "$ori: total reads input:$count\tlength good:$good\twith wrong length reads:$bad_length\n";
	
	return \%prot_cluster;

}

#-----------


sub TranslateDNASeq
   {
   	use Bio::Seq;
	my $dna=shift;
	my $l_dna = length $dna;
	# to make sequence triple
	if ($l_dna % 3  != 0){
		substr($dna,-($l_dna % 3)) = '';
	}
	my $seqobj=Bio::Seq->new(-seq =>$dna, -alphabet =>'dna');
	my $prot_obj = $seqobj->translate( 
		                                 -codontable_id => $codon,
		                                 -terminator => 'U',
		                                 -unknown => '_',
		                                 -frame => 0
		                                 );
	my $pep = $prot_obj -> seq();

}

sub COI_check{
	my $contig = shift;
	my $for_trim = $trim + 25 + 1;
	my $rev_trim = $trim + 26;
	substr($contig,0,$for_trim) = '';
	substr($contig,-$rev_trim) = '';
	(TranslateDNASeq($contig) =~ /U/  || TranslateDNASeq($contig) =~ /_/) ? return 1 : return 0;
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

			delete $depth{'N'} if (exists $depth{'N'});

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
	substr($s1,0,$l) = '';
	my $new = $s0.$s1;
	return $new;
}


sub mode_consensus {
	# In mode 2, this is only table.
	my @seqs = @_;
	my @matrix;
	my %cons_table;
	for my $i(0..$#seqs){
		my @a = split //,$seqs[$i];
		push @matrix,[@a];
	}
	my $consensus = &depth_table(@seqs);
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

sub report_depth {
	my ($hash,$title,$seq,$step,$ori) = @_;
	my %depth_sum;
	my @reports;
	my @four_bases = qw(A T C G);
	if ($ori eq 'f'){
		@reports = @{$$hash{$seq}};
		my $for_rest_len = $len - $step -1;
		for my $x(0..$for_rest_len){
			my $real_position = $x + 1;
			print HET "$title\t$real_position\t";
			for my $y(0..$#reports){
				$depth_sum{$reports[$y][$x]} ++;
			}

			foreach my $base (@four_bases){
				
				print HET (exists $depth_sum{$base} ?  "$base:$depth_sum{$base}\t" :  "$base:0\t");
			}
			print HET "\n";
			%depth_sum = ();
		}
		
	}else{
		@reports = @{$$hash{$seq}};
		my $rev_rest_len = $len -1;
		for my $x($step..$rev_rest_len){
			my $real_position = $len - $step + $x + 1;
			print HET "$title\t$real_position\t";
			for my $y(0..$#reports){
				$depth_sum{$reports[$y][$x]} ++;
			}

			foreach my $base (@four_bases){
				
				print HET (exists $depth_sum{$base} ?  "$base:$depth_sum{$base}\t" :  "$base:0\t");
			}
			print HET "\n";
			%depth_sum = ();
		}

	}
	
}
 
sub remove_outliner {
	my $hash = shift;
	my @sorted_prot_clu = sort {$#{$$hash{$b}} <=> $#{$$hash{$a}}} keys %$hash;
	for my $i(0..$#sorted_prot_clu){
		#my $size = $#{$$hash{$sorted_prot_clu[$i]}} + 1;
		my $match_id = &match($sorted_prot_clu[0],$sorted_prot_clu[$i]);
		#print "$sorted_prot_clu[$i]\t$match_id\t$size\n";
		# if in "protein" level they have more than 75% similarity, merge them together.
		if ($match_id >= $cid)
		{
			push @{$$hash{$sorted_prot_clu[0]}}, @{$$hash{$sorted_prot_clu[$i]}} ;
			# release memory
			delete $$hash{$sorted_prot_clu[$i]};
		}
	}
	return \@{$$hash{$sorted_prot_clu[0]}};
}
