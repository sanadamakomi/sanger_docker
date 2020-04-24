# sanger_docker
A shiny app packed with docker to analyse sanger sequencing data. 

Input position and ab1 file, output variants and screenshot around position.

## Before Install

### download blastn
A software named 'blastn' need be added to directory '/app', you can download the lastest version from ftp

`wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz`

### get reference

A reference file can be downloaded

`wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz`

### create index file

Three index file 'human_g1k_v37_decoy.nhr' 'human_g1k_v37_decoy.nin' & 'human_g1k_v37_decoy.nsq' can be created by 
