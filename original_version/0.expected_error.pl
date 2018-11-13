
#!/usr/bin/perl -w
use strict;
unless (@ARGV >0) {
	print "usage:perl $0 <*.fastq>";
	exit;
}
my $filename = shift;

	if ($filename =~ /.gz$/) 
{
	open FQ,"gzip -dc $filename|" or die "Can not open $filename:$!\n";
}else{
	open FQ,"$filename" or die "Can not open $filename:$!\n";
}

#Read sequences.

while(my $t1=<FQ>) {
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

	my ($aver_q,$aver_e) = &cal_ascii($qual);
	print "$t1\t$aver_q\t$aver_e\n";
}
close FQ;

#Calculate phred score.
sub cal_ascii{
	my $ascii_str = shift ;
	my $len = length $ascii_str;
	my ($exp_e,$phred) = (0,0);
	chomp $ascii_str;
	my @tmp = split //,$ascii_str;
	my @all_ascii = map {ord($_) - 33} @tmp;
	for my $i(0..$#all_ascii){
		$exp_e += 10 **(-$all_ascii[$i]/10);
		$phred += $all_ascii[$i];

	}
	#my $average_e = sprintf("%.2f",$exp_e/$len);
	my $average_q = sprintf("%.2f",$phred/$len);

	return ($average_q,$exp_e);
}

