import gzip
import subprocess


priwc = system(less -S $Pri|wc -l`)
die "the primer file need to have each forward and reverse primer"
unless ($priwc%4 == 0);

plenf = ''
plenr = ''
pris = {}
primerF = ''
primerR = ''
indp = {}
FH = {}

with open(filename,'rt') as p:
    for i in p.readlines():
        i = i.strip()
        arr = i.split()
        if len(arr) != 2:
            print("primer set is not well-formated")
        ipr = arr[1]
        if arr[0] in indp.keys():
            print("double exits $a[0] in primer set")
        else:
            ori = a[0][0,3]
            num = a[0][-3]
            $pris{$num}{$ori} = $ipr;
            pris[num][ori] = ipr
            indp[$a[0]] = ipr
            FH[a[0]] = a[0]

            if (a[0] =~ /for/i){
                $plenf = length $ipr;
                $primerF = $ipr; substr($primerF,0,$indexl) = "";

            }

            if (a[0] =~ /rev/i){
                $plenr = length $ipr;
                $primerR = $ipr; substr($primerR,0,$indexl) = "";
        
            }