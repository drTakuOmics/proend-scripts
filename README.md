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

Download (Reviewed Swiss-Pro) Arabidopsis [proteome from
Uniprot](https://www.uniprot.org/uniprotkb?query=arabidopsis&facets=reviewed%3Atrue)

1.  Fasta linearization

``` bash
zcat Arabidopsis_uniprot_proteome/uniprotkb_arabidopsis_AND_reviewed_true_2024_10_08.fasta.gz | head -n 2
```

    ## >sp|A0A067YMX8|XTH8_DIOKA Xyloglucan endotransglucosylase protein 8 OS=Diospyros kaki OX=35925 GN=XTH8 PE=1 SV=1
    ## MAASPYSIFAVQLLLLASWMLSSSSSNFNQDFNIAWGGGRARILNNGELVTLSLDKASGS

``` bash
bash HbYX-First.awk Arabidopsis_uniprot_proteome/uniprotkb_arabidopsis_AND_reviewed_true_2024_10_08.fasta.gz Arabidopsis_uniprot_proteome/arabidopsis_uniprot_proteome.ol.fa
head -n 2 Arabidopsis_uniprot_proteome/arabidopsis_uniprot_proteome.ol.fa
```

    ## >sp|A0A067YMX8|XTH8_DIOKA Xyloglucan endotransglucosylase protein 8 OS=Diospyros kaki OX=35925 GN=XTH8 PE=1 SV=1
    ## MAASPYSIFAVQLLLLASWMLSSSSSNFNQDFNIAWGGGRARILNNGELVTLSLDKASGSGFRSKNLYLFGKIDMQLKLVPGNSAGTVTTYYLSSEGSVRDEIDFEFLGNLTGEPYTLHTNVYSHGKGEREQQFRLWFDPAADFHTYSILWNSKTIVFYVDQTPVREFKNMESIGVPYLRQPMRLFSSIWNADEWATRGGLIKTDWTQAPFTTSYRNFRADNACVWAAKASSCGLAAGGNAWLSVELDAKSRGRLRWVRRNQMIYDYCVDGKRFPRGVPPECKLNLHI

2.  HbYX motif prediction

``` bash
bash HbYX-Second-Script.awk  Arabidopsis_uniprot_proteome/arabidopsis_uniprot_proteome.ol.fa  Arabidopsis_uniprot_proteome/arabidopsis_HbyX_proteome.txt
```

Total number of HbYX motif candidates:

``` bash
wc -l Arabidopsis_uniprot_proteome/arabidopsis_HbyX_proteome.txt 
```

    ## 480 Arabidopsis_uniprot_proteome/arabidopsis_HbyX_proteome.txt

3.  Testing CDC48A as HbYX proteasome regulatory protein

<!-- CDC48A (At3g09840) from Arabidopsis -->

``` bash
grep -i cdc48 -A 1 --no-group-separator Arabidopsis_uniprot_proteome/arabidopsis_HbyX_proteome.txt 
```

    ## >sp|P54609|CD48A_ARATH Cell division control protein 48 homolog A OS=Arabidopsis thaliana OX=3702 GN=CDC48A PE=1 SV=1
    ## LYN
    ## >sp|Q9SCN8|CD48D_ARATH Cell division control protein 48 homolog D OS=Arabidopsis thaliana OX=3702 GN=CDC48D PE=1 SV=1
    ## LYS
    ## >sp|Q9LZF6|CD48E_ARATH Cell division control protein 48 homolog E OS=Arabidopsis thaliana OX=3702 GN=CDC48E PE=2 SV=2
    ## LYS

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
