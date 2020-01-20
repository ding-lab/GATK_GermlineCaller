# https://hub.docker.com/r/broadinstitute/gatk
FROM broadinstitute/gatk
MAINTAINER Matthew Wyczalkowski <m.wyczalkowski@wustl.edu>

WORKDIR /opt

RUN apt-get update && apt-get install -y \
    git \
    parallel \
    vim-tiny

COPY . /opt/GATK_GermlineCaller
