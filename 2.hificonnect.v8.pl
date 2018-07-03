#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);


=head1 Name

	hificonnect.pl

=head1 Description
	
	Due to connect HIFI barcode sequence by overlaping two consensus sequences which generated 
	from Kmer method.  You can define the length of overlap, kmer length, and how many reads 
	used to make clusters, and whether to check codon translation for PCG.

=head1 Version
	
	 Version: 0.8 ( Jule 1, 2018)
	 Author: yangchentao, yangchentao@genomics.cn

=head1 Usage 
	
	perl hifi.se400.pl [options] -list fq.list 

	--list     <str>   input file, fastq file list.
	--k        <int>   kmer length,less than 21, [4];
	--min      <int>   minimun length of overlap [60]
	--max      <int>   maximum length of verlap [90]
	--oid      <float> cutoff of identity of overlap region [0.90]
	--limit    <int>   whether to limit sequencs number to save memory. 
	--index|i  <int>   index sequence lenght [4]
	--len      <int>   standard reads length [400]
	-out|o     <str>   output to file name [contig.fa]
	--cc       <?>    whether to check final COI contig's condon translation [1]
	--condon   <int>   codon table using to check translation [5], by the way, table 4,5 have same effect for COI gene.
	--frame    <0|1|2> translation start shift [1]

=head1 Change logs

	2. using own scripts instead of vsearch program to cluster reads ;

=cut

my ($help,$list,$kmer,$min,$max,$reads_lim,$overlap_id,$len,$trim,$codon,$frame,$output,$coi_check);

GetOptions(
	'help|h!' => \$help,
	'list=s' => \$list,
	'kmer|k=i' => \$kmer,
	'min:i'   => \$min,
	'max:i'   => \$max,
	'oid:f'    => \$overlap_id,
	'limit|lim:s' => \$reads_lim,
	'len:i' => \$len,
	'index|i:i' => \$trim,
	'cc!' => \$coi_check,
	'codon:i' => \$codon,
	'frame:i' => \$frame,
	'out|o:s'  => \$output
);


if (!$list or !$kmer or $help ){ 
	die `pod2text $0`;
}
if ($kmer > 21){
	die "kmer must be less than 21\n";
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


my $kmer_num = $len - $kmer + 1;
my $kmer_indx = $len - $kmer;


if ($coi_check) {
	open COI, ">$output.checked";
}
print STDERR "$output exists! overwriting ...\n" if (-e $output);

my $list_format_info = "your list is not well formated!
for example:

/path/test_For001.fastq\t/path/test_Rev001.fastq";

#`nohup sh $Bin/cpu_mem.record.sh $$ >> cpu_mem.stat$$ &`;

open OUT, ">$output";
open LOG, ">$output.log";
open HET, ">$output.depth";
open LI,"$list" or die "No such file $list!";

print LOG "## kmer = 4\n";
print LOG "## overlaping identity = $overlap_id\n";
print LOG "## overlaping step, min = $min; max = $max\n";

if ($reads_lim){
	print LOG "## reads input limitation: $reads_lim\n";
}else{
	print LOG "## Using all reads to make consensus, no limitation\n";
}

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
	my $nammme = substr($outname,-3);
	print LOG "\n//processing $outname\n";

	my ($kmers_for,$kdep_f) = &read_file($for,'f');
	print LOG ">$nammme\_f\n$kmers_for\n";
	# reverse reads have been reversed and complementioned
	my ($kmers_rev,$kdep_r) = &read_file($rev,'r');
	$kmers_rev = &comrev($kmers_rev);
	print LOG ">$nammme\_r\n$kmers_rev\n";


	# here $table_f and $table_r are two quotes of two HASH tables, key is consensus sequence
	# value is a 2D array, including each cluster's bases.

	
	#-----------anchoring overlap site--------	

	my $read0 = substr($kmers_for,-$max);
	my $read1 = substr($kmers_rev,0,$max);

	# finding overlaps position
	my %overlaps;

TNT:for (my $s=$min; $s <= $max; $s++){
		my $l0 = substr($read0,-$s);
		my $l1 = substr($read1,0,$s);
		my $id = &match($l0,$l1);

		if ($id == 1) {

			$overlaps{$s} = 1;
			print LOG "---> find position $s, with 100% mathch, jumping out loop!\n";
			last TNT; # find best result, so exit loop

		}elsif($id >= $overlap_id ){
			
			$overlaps{$s} = $id ;
		}
		
	}

	# find best overlaping result in all potenial positions
	my @candidates = sort {$overlaps{$b} <=> $overlaps{$a}} keys %overlaps;

	if (@candidates == 0 ){
		print LOG "<< no result!\n" ;
	}else  {
		my $potenial = $candidates[0];
		my $correct;
		my $makeup_consensus;
		my $s0 = substr($read0,-$potenial);
		my $s1 = substr($read1,0,$potenial);

		my @res0 = split //, $s0;
		my @res1 = split //, $s1;
		my $count = 0;
		
		#$tmp_loca0 means 0-based position in reads
		#$tmp_loca0 means 1-based position in reads

		print HET "//$nammme\t---overlap start---\n";

		for my $p(0..$#res0){
			my $tmp_loca0 = ($len - $potenial) + $p ;
			my $tmp_loca1 = $tmp_loca0 + 1;
	
			my $info;
			print HET "$nammme\t$tmp_loca1\t";
			if ($res0[$p] eq $res1[$p]){
				
				$correct .= $res0[$p];	
				$info = "$res0[$p] == $res1[$p]";
			
			}else{
				
				my $for_depth = $$kdep_f[$tmp_loca0-3] + $$kdep_f[$tmp_loca0-2] + $$kdep_f[$tmp_loca0-1] + $$kdep_f[$tmp_loca0] ;
				my $rev_depth = $$kdep_r[$p-3] + $$kdep_r[$p-2] + $$kdep_r[$p-1] + $$kdep_r[$p] ;
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
		print HET "//$nammme\t---overlap start---\n";

		$makeup_consensus = substr($kmers_for,0,$len-$potenial) . $correct . substr($kmers_rev,$potenial-$len);
		my $len_makeup_consensus = length $makeup_consensus;
		my $this_oid = $overlaps{$potenial} * 100;
		$this_oid = sprintf("%.2f",$this_oid);
		if ($coi_check){
			if (COI_check($makeup_consensus)){
				print OUT ">$nammme;overPos=$potenial;oid=$this_oid%;len=$len_makeup_consensus\n$makeup_consensus\n";
				print COI ">$nammme;overPos=$potenial;oid=$this_oid%;len=$len_makeup_consensus;check=TRUE\n$makeup_consensus\n";
			}
		}else{
			print OUT ">$nammme;overPos=$potenial;oid=$this_oid%;len=$len_makeup_consensus\n$makeup_consensus\n";
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
	my @kmers=[];
	
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

		for my $pp(0..$kmer_indx){
			push @{$kmers[$pp]}, substr($seq,$pp,$kmer);
		}
			
		last if ( $reads_lim && $count >= $reads_lim);
	}
	close FQ;

	print LOG "$ori: total reads input:$count\tlength good:$good\twith wrong length reads:$bad_length\n";

	print LOG "$ori graph:\t";
	#kmer abundance
	my $graph;
	my $path;
	my $abun_path;
	my $info;
	my @kmer_depth;
	for my $k(0..$kmer_indx){
		my %kmer_type;
		map {$kmer_type{$_}++} @{$kmers[$k]};
		my @sorted = sort{$kmer_type{$b} <=> $kmer_type{$a} or $a cmp $b} keys %kmer_type;
		if ($k == 0) {
			$graph = $sorted[0];
			$kmer_depth[0] = $kmer_type{$sorted[0]};
			$path .= "1";
			$abun_path .= "$kmer_type{$sorted[0]}";
		}else{
			STEP:for my $m(0..$#sorted){
				my $clade = $m + 1;
				if (substr($graph,(1-$kmer)) eq substr($sorted[$m],0,$kmer-1)){
					$graph .= substr($sorted[$m],-1);
					$path .= ",$clade";
					$abun_path .=",$kmer_type{$sorted[$m]}";
					$kmer_depth[$k] = $kmer_type{$sorted[$m]};
					$info .= "$k has low depth;" if ($kmer_type{$sorted[$m]} < $kmer_type{$sorted[0]}/10);
					last STEP;
				}
			}

		}
	}
	print LOG "$path\n";
	print LOG "$abun_path\n";
	if ($info){
		print LOG "[warning]: $info\n";
	}else{
		print LOG "[warning]: nono\n";
	}

	return ($graph,\@kmer_depth);

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

