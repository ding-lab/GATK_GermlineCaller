# https://hub.docker.com/r/broadinstitute/gatk
FROM broadinstitute/gatk
MAINTAINER Matthew Wyczalkowski <m.wyczalkowski@wustl.edu>

WORKDIR /opt

RUN apt-get update && apt-get install -y \
    git \
    parallel \
    vim-tiny

RUN conda install -c bioconda bcftools
RUN printf 'will cite\n' | parallel --citation

COPY . /opt/GATK_GermlineCaller
