#!/bin/bash

awk '/^>/{printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' UP000002199_224325.fasta | grep -v ^$ > multifasta.ol.fa
# Explanation:
# awk '/^>/ {...}': Process the input FASTA file
#   '/^>/': If line starts with ">", it's a header
#   'printf("\n%s\n",$0)': Print newline, then header, then another newline
#   'next': Skip to next line
#   '{ printf("%s",$0);}': For sequence lines, print without newline
#   'END {printf("\n");}': Print final newline after processing all lines
# 'UP000002199_224325.fasta': Input file name
# '| grep -v ^$': Pipe output to grep, removing empty lines
# '> multifasta.ol.fa': Redirect final output to this file
