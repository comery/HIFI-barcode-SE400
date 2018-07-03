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

	reads -> store in HASH(read=>abundance,read=>average of expected error) -> find each possible overlap ->output

=head1 Version
	
	 Version: 0.7 ( Jule 20, 2018)
	 Author: yangchentao, yangchentao@genomics.cn

=head1 Usage 
	
	perl hifi.se400.pl [options] -list fq.list 

	--list     <str>   input file, fastq file list.
	--min      <int>   minimun length of overlap [60]
	--max      <int>   maximum length of verlap [90]
	--oid      <float> cutoff of identity of overlap region [0.85]
	--cid      <float> clustering identity rate [0.95]
	--tp       <int>   how many clusters using in assembly. [2]
	--limit    <y|n>   whether to limit sequencs number to save memory. [y|n]
	--seqs     <int>   reads number limitation. [50000]
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

my ($help,$list,$min,$max,$limit,$reads_lim,$overlap_id,$cid,$len,$trim,$codon,$frame,$output,$coi_check);

GetOptions(
	'help|h!' => \$help,
	'list=s' => \$list,
	'min:i'   => \$min,
	'max:i'   => \$max,
	'oid:f'    => \$overlap_id,
	'cid:f'   => \$cid,
	'limit|lim:s' => \$limit,
	'seqs:i'  => \$reads_lim,
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

if ($limit eq 'y'){
	$reads_lim ||= 50000;
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

	my ($cluster_for,$abun_f) = &read_file($for,'f');

	my ($cluster_rev,$abun_r) = &read_file($rev,'r');

	
	my @consensus_for = sort {$$abun_f{$b} <=> $$abun_f{$a} or $a cmp $b} keys %$abun_f;
	my @consensus_rev = sort {$$abun_r{$b} <=> $$abun_r{$a} or $a cmp $b} keys %$abun_r;

#print "@consensus_for\n";
	print STDERR "$outname\n";
	print STDERR "Overlapping...\n";
	
	#-----------anchoring overlap site--------
	
	my $nammme = substr($outname,-3);
	my $paired = 0;
FOR:for my $i(0..$#consensus_for){
		#next if (undef $consensus_for[$i]);
		my $cluster_f = $consensus_for[$i];
		print "$consensus_for[$i]\n";
		
		my $pos1 = $i + 1;
	REV:for my $j (0..$#consensus_rev){
			next if ($consensus_rev[$j] eq '');
			my $cluster_r = $consensus_rev[$j];
			my $pos2 = $j + 1;
			my $c_size = ($$abun_f{$cluster_f} + $$abun_r{$cluster_r}) / 2;
			# print just once
			
			my $read0 = substr($cluster_f,-$max);
			my $read1 = substr($cluster_r,0,$max);

		# finding overlaps position-----------------------
			my %overlaps;

	OVERLAP:for (my $s=$min; $s <= $max; $s++){
				my $l0 = substr($read0,-$s);
				my $l1 = substr($read1,0,$s);
				my $id = &match($l0,$l1);

				if ($id == 1) {
					$paired ++;
					$overlaps{$s} = 1;
					print STDERR "$pos1-$pos2 ---> find position $s 100% mathch, jumping out loop!\n";
					delete $consensus_for[$i];
					$consensus_rev[$j] = '';
					last OVERLAP; # find best result, so exit loop

				}elsif($id >= $overlap_id ){
					
					$overlaps{$s} = $id ;
					

				}else{
					next OVERLAP;
				}
				
			}
		# overlaping finished----------------------------

			print STDERR "$pos1-$pos2 Finding overlap done\t";
			# find best overlaping result in all potenial positions
			my @candidates = sort {$overlaps{$b} <=> $overlaps{$a}} keys %overlaps;
		
			if (@candidates == 0 ){
				print STDERR "<< $pos1-$pos2 has no result!\n" ;
				next REV;

			}else  {

				print STDERR "$pos1-$pos2 connected\n";
				$paired ++;
				delete $consensus_for[$i];
				$consensus_rev[$j] = '';

				my $potenial = $candidates[0];
				my $correct;
				my $makeup_consensus;

				my $s0 = substr($read0,-$potenial);
				my $s1 = substr($read1,0,$potenial);

				my @res0 = split //, $s0;
				my @res1 = split //, $s1;
				my $count = 0;
				
				for my $p(0..$#res0){
					my $tmp_loca0 = ($len - $potenial) + $p ;
					my $for_error = ${$$cluster_for{$consensus_for[$i]}}[$tmp_loca0];
					my $rev_error = ${$$cluster_rev{$consensus_rev[$j]}}[$p];
					

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
					my $info;
					if ($res0[$p] eq $res1[$p]){
						$correct .= $res0[$p];
						$info = "$res0[$p] == $res1[$p]";
					}else{
						if ($for_error < $rev_error) {
							$correct .= $res0[$p];
							$info = "[F=$res0[$p]:$for_error < R=$res1[$p]:$rev_error]";
						}elsif($for_error == $rev_error){
							print LOG "$nammme\t$tmp_loca1\tcould be a heterozygosis site!\n";
							$correct .= $res0[$p];
							$info = "[F=$res0[$p]:$for_error = R=$res1[$p]:$rev_error]";
						}else{
							$correct .= $res1[$p];
							$info = "[F=$res0[$p]:$for_error > R=$res1[$p]:$rev_error]";
						}
					}
					
					print HET "$info\n";

				}
				
				$makeup_consensus = substr($cluster_f,0,$len-$potenial) . $correct . substr($cluster_r,$potenial-$len);
				my $len_makeup_consensus = length $makeup_consensus;
				my $this_oid = $overlaps{$potenial} * 100;
				$this_oid = sprintf("%.2f",$this_oid);
				if ($coi_check){
					if (COI_check($makeup_consensus)){
						print OUT ">$nammme\_$pos1-$pos2;$$abun_f{$cluster_f}\_$$abun_r{$cluster_r};size=$c_size;overPos=$potenial;oid=$this_oid%;len=$len_makeup_consensus\n$makeup_consensus\n";
						print COI ">$nammme\_$pos1-$pos2;$$abun_f{$cluster_f}\_$$abun_r{$cluster_r};size=$c_size;overPos=$potenial;oid=$this_oid%;len=$len_makeup_consensus;check=TRUE\n$makeup_consensus\n";
					}
				}else{
					print OUT ">$nammme\_$pos1-$pos2;$$abun_f{$cluster_f}\_$$abun_r{$cluster_r};size=$c_size;overPos=$potenial;oid=$this_oid%;len=$len_makeup_consensus\n$makeup_consensus\n";
				}

			}
			next FOR;
			
		}
	}
	my $num_f = @consensus_for;
	my $num_r = @consensus_rev;
	print LOG "forward clusters: $num_f\nreverse clusters: $num_r\n";
	print LOG "forward-reverse paired: $paired\n";
}

close OUT;
close COI;

#--------------------------------------------------------------------------

# read fastq file and remove outliers 
sub read_file{

	my $filename = shift;
	my $ori = shift;
	my %cluster;
	my %abun;
	
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
			$abun{$seq} ++;
			push @{$cluster{$seq}},$qual;
		}else{
			$abun{&comrev($seq)} ++;
			push @{$cluster{&comrev($seq)}},reverse $qual;
		}
		
			
		last if ( $limit eq 'y' && $count >= $reads_lim);
	}
	close FQ;

	print LOG "$ori: total reads input:$count\tlength good:$good\twith wrong length reads:$bad_length\n";
	

	#rebuild the hash, key is unique sequence, value is average expected 
	my @sorted_seqs = sort{$cluster{$b} <=> $cluster{$a}} keys %cluster;
	for my $i(0..$#sorted_seqs){
		if (@{$cluster{$sorted_seqs[$i]}} == 1){
			@{$cluster{$sorted_seqs[$i]}} = split //,${$cluster{$sorted_seqs[$i]}}[0];
		}else{
			@{$cluster{$sorted_seqs[$i]}} = &average_expected_e(@{$cluster{$sorted_seqs[$i]}});
		}

	}

	return (\%cluster,\%abun);

}

#-----------

#Calculate expected error of reads


sub average_expected_e
	{
		my @seqs = @_;
		my @arr_e;
		my $total_e = 0;
		my $llen = length ($seqs[0]);
		for(my $y=0;$y<$llen;$y++){
			foreach my $align(@seqs){
				my $aa=substr($align,$y,1);
				my $phred = ord($aa) - 33;
				$total_e += 10 **(-$phred/10);; 
			}
			my $average_e = sprintf("%.3f",$total_e / @seqs);
			$total_e = 0;
		
			push @arr_e, $average_e;
			
		}
		
		return @arr_e;
}

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



 
