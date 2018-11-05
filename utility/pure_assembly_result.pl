#! /usr/bin/perl -w
use strict;
open IN,shift;
my $current_sam;
my %hash;
my @all_id;
my %seqs;
my %full_name;
while (<IN>){
	chomp;
	s/^>//;
	my $info = $_;
	my @a = split /;/,$info;
	my $id = $a[0];
	$full_name{$id} = $info;
	my $size = $1 if ($info =~ /size=(.+?);/) ;
	my $sam = $1 if ($a[0] =~ /(\d+)_/);
	chomp(my $seq = <IN>);
	$seqs{$id} = $seq;
	# choose most abundant contig, and filter out short amplicons.
	if (! $current_sam) {
		$current_sam = $sam;
		$hash{$id} = $size;
	}elsif ($sam eq $current_sam) {
		$hash{$id} = $size;
	}else{
		my @sorted_key = sort{$hash{$b} <=> $hash{$a} or $a cmp $b} keys %hash;
		my ($s1,$s2) = ($1,$2) if ($full_name{$sorted_key[0]} =~ /;(\d+)_(\d+);/);
		push @all_id, $sorted_key[0] if ($s1 >= 10 && $s2 >= 10);
		%hash = ();
		$current_sam = $sam;
		$hash{$id} = $size;
	}
}
# deal last item
my @sorted_key = sort{$hash{$b} <=> $hash{$a} or $a cmp $b} keys %hash;
my ($s1,$s2) = ($1,$2) if ($full_name{$sorted_key[0]} =~ /;(\d+)_(\d+);/);
push @all_id, $sorted_key[0] if ($s1 >= 10 && $s2 >= 10);
%hash = ();

foreach my $k(@all_id){
	if (length($seqs{$k}) > 650){
		print ">$full_name{$k}\n$seqs{$k}\n";
	}
}

close IN;
