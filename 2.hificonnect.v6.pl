#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);


=head1 Name

	hificonnect.pl

=head1 Description
	
	Due to connect HIFI barcode sequence by overlaping two consensus sequences which generated 
	from clustering method, in this program I use VSEARCH.  You can define the length of overlap,
	how many reads used to make clusters, and whether to check codon translation for PCG.

	reads -> find most abundant seed and remove outliers -> call consensus -> overlapping

=head1 Version
	
	 Version: 0.6 ( Jule 13, 2018)
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

	2. using own scripts instead of vsearch program to cluster reads ;

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
$cid ||= 0.97;
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

`nohup sh $Bin/cpu_mem.record.sh $$ >> cpu_mem.stat$$ & `;

open OUT, ">$output";
open LOG, ">$output.log";
open HET, ">$output.depth";
open LI,"$list" or die "No such file $list!";

print HET "#ID\tPos\tA\tT\tC\tG\tfor vs rev\n";
print LOG "##overlaping identity = $overlap_id\n";
print LOG "##clustering similarity = $cid\n";
print LOG "##overlaping step, min = $min; max = $max\n";

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
	print LOG "\n//processing $outname\n";

	my $cluster_for = &read_file($for,'f');
# reverse reads have been reversed and complementioned
	my $cluster_rev = &read_file($rev,'r');


	# here $table_f and $table_r are two quotes of two HASH tables, key is consensus sequence
	# value is a 2D array, including each cluster's bases.
	
	my ($table_f,$consensus_f,$abundance_f) = &make_depth_table($cluster_for);

	my ($table_r,$consensus_r,$abundance_r) = &make_depth_table($cluster_rev);

	my $c_size = ($abundance_f + $abundance_r) / 2;

	my @consensus_for = ($consensus_f);
	my @consensus_rev = ($consensus_r);

	print STDERR "$outname\n";
	print STDERR "Making consensus done\n";
	print STDERR "Overlapping...\n";
	
	#-----------anchoring overlap site--------
	
	my $nammme = substr($outname,-3);

	for my $i(0..$#consensus_for){
		my $cluster_f = $consensus_for[$i];
		my $pos1 = $i + 1;
		print LOG ">$nammme\_f\_$pos1\tsize=$abundance_f\n$cluster_f\n";
		for my $j (0..$#consensus_rev){
			my $cluster_r = $consensus_rev[$j];
			my $pos2 = $j + 1;
			# print just once
			print LOG">$nammme\_r\_$pos2\tsize=$abundance_r\n$cluster_r\n" if ($pos1 == 1);
			
			my $read0 = substr($cluster_f,-$max);
			my $read1 = substr($cluster_r,0,$max);

		# finding overlaps position
			my %overlaps;

		TNT:for (my $s=$min; $s <= $max; $s++){
				my $l0 = substr($read0,-$s);
				my $l1 = substr($read1,0,$s);
				my $id = &match($l0,$l1);

				if ($id == 1) {

					$overlaps{$s} = 1;
					print LOG "---> find position $s 100% mathch, jumping out loop!\n";
					last TNT; # find best result, so exit loop

				}elsif($id >= $overlap_id ){
					
					$overlaps{$s} = $id ;
				}
				
			}

			print STDERR "$pos1-$pos2 Finding overlap done\n";
			# find best overlaping result in all potenial positions
			my @candidates = sort {$overlaps{$b} <=> $overlaps{$a}} keys %overlaps;
			
		
			if (@candidates == 0 ){
				print LOG "<< $pos1-$pos2 has no result!\n" ;
			}else  {
				my $potenial = $candidates[0];
				my $correct;
				my $makeup_consensus;

				my $s0 = substr($read0,-$potenial);
				my $s1 = substr($read1,0,$potenial);

				#report forward depth
				&report_depth($table_f,"$nammme\_$pos1-$pos2\t",$potenial,'f');
				#report depth of overlap
				print HET "$nammme\_$pos1-$pos2\t\t-----overlap start-----\n";

				my @res0 = split //, $s0;
				my @res1 = split //, $s1;
				my $count = 0;
				
				for my $p(0..$#res0){
					my $tmp_loca0 = ($len - $potenial) + $p ;
					my $for_depth = ${$$table_f[$tmp_loca0]}{$res0[$p]};
					my $rev_depth = ${$$table_r[$p]}{$res1[$p]};
					

					# site is changed, be careful!
					# !!!!
						

#forward read
#								$q									
#0	1	2	3	4	5	6	7	8	9	10									forward length = 11
#A	T	C	G	G	A	A	G	T	T	A									reverse length =11
#							G	T	A	A	A	C	G	C	G	G	A		overlap = 4
#							0	1	2	3	4	5	6	7	8	9	10		$q = for_len - verlap + $p
#								$p											8 = 11 - 4 + 1
#reverse read

					my $tmp_loca1 = $tmp_loca0 + 1;
					print HET "$nammme\_$pos1-$pos2\t$tmp_loca1\t";
					my %sum;
					my @four_bases = qw(A T C G);
					foreach (@four_bases){
						exists ${$$table_f[$tmp_loca0]}{$_} ? 1 : (${$$table_f[$tmp_loca0]}{$_} = 0);
						exists ${$$table_r[$p]}{$_} ? 1 : (${$$table_r[$p]}{$_} = 0);
						$sum{$_} = ${$$table_f[$tmp_loca0]}{$_} + ${$$table_r[$p]}{$_} ;
						print HET "$_:$sum{$_}\t";
					}
					#map {$sum{$_} = ${$$table_f[$tmp_loca0]}{$_} + ${$$table_r[$p]}{$_}} @four_bases;

					%sum = ();
					
					my $info;
					if ($res0[$p] eq $res1[$p]){
						$correct .= $res0[$p];
						$info = "$res0[$p] == $res1[$p]";
					}else{
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
					}
					
					print HET "$info\n";
				}

				print HET "$nammme\_$pos1-$pos2\t\t-----overlap end-----\n";
				
				#report reverse depth
				&report_depth($table_r,"$nammme\_$pos1-$pos2\t",$potenial,'r');
				
				$makeup_consensus = substr($cluster_f,0,$len-$potenial) . $correct . substr($cluster_r,$potenial-$len);
				my $len_makeup_consensus = length $makeup_consensus;
				my $this_oid = $overlaps{$potenial} * 100;
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

#--------------------------------------------------------------------------

# read fastq file and remove outliers 
sub read_file{

	my $filename = shift;
	my $ori = shift;
	my %cluster;
	
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

		if ($ori eq 'f'){
			$cluster{$seq} ++;
		}else{
			$cluster{&comrev($seq)} ++;
		}
		
			
		last if ( $reads_lim && $count >= $reads_lim);
	}
	close FQ;

	print LOG "$ori: total reads input:$count\tlength good:$good\twith wrong length reads:$bad_length\n";
	

	# remove outlier

	my @sorted_seqs = sort{$cluster{$b} <=> $cluster{$a}} keys %cluster;
	#print LOG "top1 read has low abundance!\n" if ($cluster{$sorted_seqs[0]} < 10);
	print LOG ">top1_seed\tsize=$cluster{$sorted_seqs[0]}\n$sorted_seqs[0]\n";
	for my $ss (1..$#sorted_seqs){
		my $simi = &match($sorted_seqs[0],$sorted_seqs[$ss]);
		delete $cluster{$sorted_seqs[$ss]} if ($simi < $cid) ;
	}

	return \%cluster;

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
		my $hash = shift;
		my @seqs = keys %$hash;
		my %depth;
		my $consensus;
		my $llen = length ($seqs[0]);
		for(my $y=0;$y<$llen;$y++){
			foreach my $align(@seqs){
				my $aa=substr($align,$y,1);
				$depth{$aa} += $$hash{$align}; ## $a = (A|T|C|G|N);
			}

			delete $depth{'N'} if (exists $depth{'N'});

			my @sort_base = sort{$depth{$b} <=> $depth{$a} or $a cmp $b} keys %depth;
		
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


sub make_depth_table {

	#summary quality and depth of each kind base
	my @base_sum;
	my $consensus;
	my $size = 0;
	my $hash = shift;
	my @seqs = keys %$hash;
	
	for my $i(0..$len -1){ 
		my %bases_distribution;
		foreach my $j(@seqs){
			my $tmp = substr($j,$i,1);
			$bases_distribution{$tmp} += $$hash{$j};
			$size += $$hash{$j} if ($i == 0);
		}

		delete $bases_distribution{'N'} if (exists $bases_distribution{'N'});

		my @sort_base = sort{$bases_distribution{$b} <=> $bases_distribution{$a} or $a cmp $b} keys %bases_distribution;

		# for some case, different bases have same deepth
	
		$consensus .= $sort_base[0];

		$base_sum[$i] = \%bases_distribution;
	}

	return (\@base_sum,$consensus,$size);
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
	my ($arr,$title,$step,$ori) = @_;
	my @four_bases = qw(A T C G);

	if ($ori eq 'f'){
		my $for_rest_len = $len - $step -1;
		for my $x(0..$for_rest_len){
			my $real_position = $x + 1;
			print HET "$title\t$real_position\t";
			foreach my $base (@four_bases){
				
				print HET (exists ${$$arr[$x]}{$base} ?  "$base:${$$arr[$x]}{$base}\t" :  "$base:0\t");
			}
			print HET "\n";
		
		}
		
	}else{
		my $rev_rest_len = $len -1;
		for my $x($step..$rev_rest_len){
			my $real_position = $len - $step + $x + 1;
			print HET "$title\t$real_position\t";

			foreach my $base (@four_bases){
				
				print HET (exists ${$$arr[$x]}{$base} ?  "$base:${$$arr[$x]}{$base}\t" :  "$base:0\t");
			}
			print HET "\n";
		}

	}
	
}
 
