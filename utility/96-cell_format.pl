# make 96-cell plate format
my @b = [];
open IN,shift;
while(<IN>){
	chomp;
	my @aa = split;
	my $id = substr($aa[0],-3);
	my $int = int(($id-1)/8);
	my $res = ($id-1)%8 + 1;
	$b[$res][$int] = "$aa[1]\t$aa[3]";

}


for my $x (0 .. $#b){
	print join ("\t",@{$b[$x]});
	print "\n";
    
}
