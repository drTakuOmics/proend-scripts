#!/bin/bash

awk ' BEGIN {OFS = "\n"} {header = $0 ; getline seq ; if (seq ~ /[G|L|I|V|A|Y|F|W|C|M ]Y[A-Z]$/){ taa=substr(seq,length(seq)-2,length(seq)) ;  print header"\n"taa}}' multifasta.ol.fa > 2aa.txt
# Explanation:
# 'BEGIN {OFS = "\n"}': Set output field separator to newline
# '{header = $0 ; getline seq': Store current line as header, read next line as sequence
# 'if (seq ~ /[G|L|I|V|A|Y|F|W|C|M ]Y[A-Z]$/)': Check if sequence ends with HbYX-motif
# 'taa=substr(seq,length(seq)-2,length(seq))': Extract last 3 amino acids
# 'print header"\n"taa': Print header and extracted motif
# 'multifasta.ol.fa': Input file (replace with your FASTA file if already in one-line-per-sequence format)
# '> 2aa.txt': Output results to this file

# Note: If your FASTA file is already in one-line-per-sequence format,
# simply replace 'multifasta.ol.fa' with the name of your FASTA file and execute this script.
