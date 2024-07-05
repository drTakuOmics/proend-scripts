#!bin/bash/ 
awk ' BEGIN {OFS = "\n"} {header = $0 ; getline seq ; if (seq ~ /[G|L|I|V|A|Y|F|W|C|M ]Y[A-Z]$/){ taa=substr(seq,length(seq)-2,length(seq)) ;  print header"\n"taa}}'  multifasta.ol.fa  > 2aa.txt
