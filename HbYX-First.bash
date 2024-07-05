#!/bin/bash

awk '/^>/{printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' UP000002199_224325.fasta | grep -v ^$ > multifasta.ol.fa

