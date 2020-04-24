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

# 中文版说明

本项目基于shiny app，能可视化地分析一代序列，输入位点和ab1文件，可以看色谱图，call变异并自动截图，为了方便部署使用docker打包shiny app。

## 安装前

除了需要git clone本项目之外，还需要下载blastn和参考基因组文件，构建index。因为这些文件太大了，就没有写在Dockerfile里。

### Blastn

可以去ncbi的ftp获得最新的blast，本项目在构建时使用的是2.9.0版：

`wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz`

解压后进入bin目录，拷贝blastn文件到本项目的app目录下。

### 参考基因组

下载GRch37参考基因组：

`wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz`

然后构建比对需要的索引文件：

`./ncbi-blast-2.9.0+/bin/makeblastdb -in human_g1k_v37_decoy.fasta -dbtype nucl -parse_seqids -out ./human_g1k_v37_decoy`

然后拷贝'human_g1k_v37_decoy.nhr', 'human_g1k_v37_decoy.nin' & 'human_g1k_v37_decoy.nsq'三个文件到app目录下。

## 创建镜像

`docker build -t sanger_docker .`

## 运行容器

可以指定一个端口，比如80

`docker run -p 80:3838 sanger_docker`

