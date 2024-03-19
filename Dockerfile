# centos7
FROM rstudio/r-base:4.3-centos7

# system libraries of general use
RUN yum -y install epel-release && \
    yum update -y && yum install -y \
    sudo \
    pandoc \
    poppler-cpp-devel \
    ImageMagick-c++-devel \
    qpdf

# basic shiny functionality
RUN R -e "install.packages(c('shiny', 'rmarkdown'), repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"

# install dependencies of the omim app
RUN R -e "install.packages(c('shinythemes', 'shinydashboard', 'knitr', 'stringr', 'pdftools', 'magick', 'digest', 'DT', 'BiocManager', 'RCurl'), repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"
RUN R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');BiocManager::install(c('S4Vectors', 'XVector', 'Biostrings', 'GenomeInfoDb', 'IRanges', 'GenomicRanges', 'GenomeInfoDb', 'sangerseqR'))"

# copy the app to the image
RUN mkdir -p /var/data && \
    mkdir -p /var/data/log

COPY app /var/data
COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

WORKDIR /var/data

CMD ["R", "-e", "shiny::runApp('/var/data')"]
