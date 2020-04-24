# Install R version 3.6
FROM rocker/r-ver:3.6.0

# Use tsinghua mirrors
RUN mv /etc/apt/sources.list /etc/apt/sources.list.bak && \
    echo "deb http://mirrors.ustc.edu.cn/debian/ stretch main contrib non-free" >> /etc/apt/sources.list && \
    echo "deb http://mirrors.ustc.edu.cn/debian/ stretch-updates main contrib non-free" >> /etc/apt/sources.list && \
    echo "deb http://mirrors.ustc.edu.cn/debian/ stretch-backports main contrib non-free" >> /etc/apt/sources.list && \
    echo "deb http://mirrors.ustc.edu.cn/debian-security/ stretch/updates main contrib non-free" >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list


# Install Ubuntu packages
# libmagick++-dev & libpoppler-cpp-dev are required for pdftools & magick
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget \
    libpoppler-cpp-dev \
    libmagick++-dev
            

# Install shiny
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('shiny', 'shinythemes', 'shinydashboard', 'knitr', 'rmarkdown', 'stringr', 'pdftools', 'magick', 'digest', 'DT', 'BiocManager'), repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')" && \
    R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');BiocManager::install(c('IRanges', 'GenomicRanges', 'GenomeInfoDb', 'sangerseqR'))"


# Copy configuration files into the Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY /app /srv/shiny-server/


# Make the ShinyApp available at port 3838
EXPOSE 3838


RUN chmod +x /usr/bin/shiny-server.sh && \
    chmod -R 777 /srv/shiny-server && \
    mkdir -p /project && \
    chown -R shiny /project

CMD ["/usr/bin/shiny-server.sh"]
