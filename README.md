ProEnd Pipeline
================
David & Aimer

# ProEnd Scripts

This project contains two Bash scripts designed to handle and analyze
multiple protein sequences. The scripts streamline the extraction and
identification of specific C-terminal protein HbYX-motifs from a given
FASTA file. The third script shows how to download the proteomes from
UniProt.

If you have just one protein sequence or your file is already one
protein sequence per line go to script 2 \## Scripts

### 1. HbYX-First.bash

This script formats a multi-entry FASTA file into a single line per
entry format, preparing it for further analysis. Use this script first
if you have multiple sequences, a proteome, or alignments.

#### Usage

Ensure you have the necessary FASTA file in the same directory or
specify the path to the file. Execute the script by running:

`#!bash ./HbYX-First.bash`

or

`bash HbYX-First.bash` This will output a file named `multifasta.ol.fa`,
containing all the sequences from the original FASTA file, formatted for
further processing.

### 2. HbYX-Second-Script.bash

This script searches for a specific HbYX-motif at the C-terminus of the
protein sequences in the `multifasta.ol.fa` file created by the first
script.

#### Usage

Run the script using:

`#!bash ./HbYX-Second-Script.bash`

or

`bash HbYX-Second-Script.bash` It will produce a file named `2aa.txt`,
containing the motifs found along with the corresponding header from the
FASTA file, if available.

### Testing ProEND with *Arabidopsis proteome*

## Requirements

- Unix-like environment
- AWK installed

## Installation

No installation is required. Simply clone this repository or download
the scripts to your local machine.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Data Folder HbYX_data_tables

This folder contains data results files for the ProEnd Scripts project.

## Cite

This code can be cited currently as BioRxiv preprint
\[[1](#ref-salcedo2024proend)\]

## References

<div id="refs" class="references csl-bib-body">

<div id="ref-salcedo2024proend" class="csl-entry">

1\. Salcedo-Tacuma DM, Howells G, McHose C, Gutierrez-Diaz A, Smith DM.
ProEnd: A comprehensive database for identifying HbYX motif-containing
proteins across the tree of life. bioRxiv. 2024;2024–06.

</div>

</div>
