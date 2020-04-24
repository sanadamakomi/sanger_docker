# sanger_docker
A shiny app packed with docker to analyse sanger sequencing data. 

Input position and ab1 file, output variants and screenshot around position.

## Before Install

### Blastn

You can download the lastest version from ftp. This app use version 2.9.0:

`wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz`

Depress and copy `blastn` to the directory `/app`.

### Database for blast

A reference file can be downloaded from ftp:

`wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz`

Three index file  can be created by:

`./ncbi-blast-2.9.0+/bin/makeblastdb -in human_g1k_v37_decoy.fasta -dbtype nucl -parse_seqids -out ./human_g1k_v37_decoy`

And copy three file 'human_g1k_v37_decoy.nhr', 'human_g1k_v37_decoy.nin' & 'human_g1k_v37_decoy.nsq' to the directory `/app`.

## Build Image

`docker build -t sanger_docker .`

## Run Container

`docker run -p 80:3838 sanger_docker`
