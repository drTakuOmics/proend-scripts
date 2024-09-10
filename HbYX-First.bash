#!/bin/bash

# This script processes a FASTA file to reformat it into a single line per sequence format

# Use awk to process the input file
awk '
    /^>/ {
        # If the line starts with ">", it's a header
        # Print a newline (except for the first header), then print the header
        printf("\n%s\n", $0);
        next;  # Skip to the next line
    }
    {
        # For sequence lines, print without newline
        printf("%s", $0);
    }
    END {
        # After processing all lines, print a final newline
        printf("\n");
    }
' UP000002199_224325.fasta |  # Input file name
    grep -v ^$ > multifasta.ol.fa  # Remove empty lines and save to output file

# Explanation of the grep command:
# -v: Invert the match, i.e., select non-matching lines
# ^$: Regex pattern for an empty line (start of line immediately followed by end of line)
# > multifasta.ol.fa: Redirect output to this file
