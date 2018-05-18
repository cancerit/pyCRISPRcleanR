FROM  ubuntu:16.04
MAINTAINER  cgphelp@sanger.ac.uk

LABEL uk.ac.sanger.cgp="Cancer Genome Project, Wellcomei Sanger Institute" \
      version="1.0.0" \
      description="Tool to perform crisprcleaner analysis"

USER root

ENV OPT /opt/wtsi-cgp
ENV PATH $OPT/bin:$PATH
ENV LD_LIBRARY_PATH $OPT/lib
ENV R_LIBS $OPT/R-lib
ENV R_LIBS_USER $R_LIBS

RUN mkdir -p $R_LIBS_USER
# install system tools
RUN apt-get update && \
  apt-get install -yq --no-install-recommends lsb-release && \
  echo "deb http://cran.rstudio.com/bin/linux/ubuntu $(lsb_release -cs)/" \
  >> /etc/apt/sources.list && \
  gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9 && \
  gpg -a --export E084DAB9 | apt-key add - && \
  apt-get update && \
  apt-get install -qy --no-install-recommends \
    libcairo2-dev \
    r-base \
    r-base-dev \
    python3  \
    python3-dev \
    python3-setuptools \
    python3-pip

RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("DNAcopy", ask=FALSE, lib="'"${R_LIBS_USER}"'")'

# install crisprcleanr
#RUN pip3 --no-cache-dir install https://github.com/cancerit/pyCRISPRcleanR/releases/download/1.1.1/pyCRISPRcleanR-1.1.1-py3-none-any.whl
COPY pyCRISPRcleanR-1.1.1-py3-none-any.whl $OPT
RUN pip3 --no-cache-dir install $OPT/pyCRISPRcleanR-1.1.1-py3-none-any.whl
### security upgrades and cleanup
RUN apt-get -yq update && \
    apt-get -yq install unattended-upgrades && \
    unattended-upgrades

RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER ubuntu

WORKDIR /home/ubuntu

CMD ["/bin/bash"]
