#! usr/bin/perl -w
use strict;
=head1 Description

        Usage: perl extract.pl <parameter>

        --fq         "xx.1.fq"
        --pri        "primer set with index ahead, format in "ID \t seq";
        --out        "output directory and prefix of the output file"
        --seqtype    "output format [fa|fq]"
        --help       "print out this information.

=cut

use Getopt::Long;
use FindBin qw($Bin $Script);
use strict;


my ($Fq,$Pri,$Out,$Help,$format,$indexl);


GetOptions(
        "fq:s"=>\$Fq,
        "out:s"=>\$Out,
        "pri:s"=>\$Pri,
		"seqtype:s"=>\$format,
         "help"=>\$Help
);
die `pod2text $0` if ($Help || !defined ($Fq) || !defined ($Pri));
$Out ||= "taged";
$format |= "fq";
$indexl ||= 4;
chomp(my $wd = `pwd`);

if ($Fq=~/gz$/){
	open (FFQ, "<:gzip",$Fq) || die $!;
}else{
    open (FFQ, $Fq) || die $!;
}

my $priwc;
$priwc = `less -S $Pri|wc -l`;
die "the primer file need to have each forward and reverse primer"
unless ($priwc%4 == 0);

#die "the primer file contain more than 2 primer sets, however, without index length hava been seted" if ($priwc > 4 );

mkdir $Out unless (-d $Out);
my (%indp,%FH);
open PRI, "$Pri" || die $!;
my ($plenf, $plenr,%pris,$primerF,$primerR);
while (<PRI>){
	chomp;
	my @a=split /\s+/;
	die "primer set is not in correct format\n" unless (@a==2);
	my $ipr = $a[-1];

	if (exists $indp{$a[0]}){
		warn "double exits $a[0] in primer set\n";
	}else{
        my $ori = substr($a[0],0,3);
        my $num = substr($a[0],-3);
        $pris{$num}{$ori} = $ipr;
		$indp{$a[0]} = $ipr;
        $FH{$ipr} = $a[0];
		if ($a[0] =~ /for/i){
			$plenf = length $ipr;
			$primerF = $ipr; substr($primerF,0,$indexl) = "";

		}

		if ($a[0] =~ /rev/i){
			$plenr = length $ipr;
			$primerR = $ipr; substr($primerR,0,$indexl) = "";
	
		}
	}
}
close PRI;
my $neg_priF = comrev($primerF);
my $neg_priR = comrev($primerR);


open (LS,">","$Out.splitlist");
my @sorted = sort {$a <=> $b} keys %pris;
for my $key (@sorted){
	print LS "$wd/$Out/$Out\_For$key\.fasta\t$wd/$Out/$Out\_Rev$key\.fasta\n" if ($format eq 'fa');
	print LS "$wd/$Out/$Out\_For$key\.fastq\t$wd/$Out/$Out\_Rev$key\.fastq\n" if ($format eq 'fq');

}
close LS;

my %filehandle;
for my $key (keys %indp){
	open ($filehandle{"$key"}, ">", "$Out/$Out\_$key\.fasta")  if ($format eq 'fa');
	open ($filehandle{"$key"}, ">", "$Out/$Out\_$key\.fastq")  if ($format eq 'fq');

}

open ERR, ">$Out\_err.fasta";


my $seqnum=0;
my $err = 0;
my $assigned = 0;
my %count;
while(my $t1=<FFQ>){
	$seqnum++;
	chomp($t1);
	$t1 =~s/^@//;
	chomp(my $seq=<FFQ>);
	<FFQ>;
	chomp(my $qual = <FFQ>);
	my $headf = substr($seq,0,$plenf);
	my $headr = substr($seq,0,$plenr);
	next if ($seq =~ /N/);
    my $tmp = $seq;
	if (exists $FH{$headf}) {
        substr($tmp,0,$plenf) = "";
		$count{$FH{$headf}}{'a'} ++;
		
		if ($tmp !~ /$primerF/ && $tmp !~ /$primerR/ && $tmp !~ /$neg_priF/ && $tmp !~ /$neg_priR/){
			if ($format eq 'fq'){print {$filehandle{$FH{$headf}}} "@";}else{print {$filehandle{$FH{$headf}}} ">";}
			print  {$filehandle{$FH{$headf}}}  "$FH{$headf}\_$seqnum\n$seq\n" ;
			print {$filehandle{$FH{$headf}}} "+\n$qual\n" if ($format eq 'fq');
			$count{$FH{$headf}}{'t'} ++;
			$assigned ++;
		}else{
			$err ++;
		}
	}elsif(exists $FH{$headr}){
        substr($tmp,0,$plenr) = "";
		$count{$FH{$headr}}{'a'} ++;
		
		if ($tmp !~ /$primerF/ && $tmp !~ /$primerR/ && $tmp !~ /$neg_priF/ && $tmp !~ /$neg_priR/){
			if ($format eq 'fq'){print {$filehandle{$FH{$headr}}} "@";}else{print {$filehandle{$FH{$headr}}} ">";}
			print  {$filehandle{$FH{$headr}}}  "$FH{$headr}\_$seqnum\n$seq\n";
			print {$filehandle{$FH{$headr}}} "+\n$qual\n" if ($format eq 'fq');
			$count{$FH{$headr}}{'t'} ++;
			$assigned ++;
		}else{
			$err ++;
		}
	}else{
		#print ERR ">$seqnum\_f\n$seq\n";
		$err ++;
	}
}	
close ERR;

open LOG,">$Out.assign.log";

print LOG "total reads:\t$seqnum\n";
print LOG "err reads:\t$err\n";
print LOG "assigned:\t$assigned\n";
foreach my $i(sort keys %count){
	print LOG "$i\t$count{$i}{'a'}\t$count{$i}{'t'}\n";
}

close LOG;
print "done\n";

sub comrev {
    my $a = shift;
    chomp $a;
    $a =~ tr/ATCG/TAGC/;
    $a= reverse $a;
    return $a;
}
