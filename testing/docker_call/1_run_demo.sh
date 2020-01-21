source ../../docker/docker_image.sh

DATD="../demo_data"
OUTD="./output"

BAM="/data/HCC1954.NORMAL.30x.compare.COST16011_region.bam"
REF="/data/Homo_sapiens_assembly19.COST16011_region.fa"

PROCESS="/opt/GATK_GermlineCaller/src/process_sample_parallel.sh"

# Using python to get absolute path of DATD.  On Linux `readlink -f` works, but on Mac this is not always available
# see https://stackoverflow.com/questions/1055671/how-can-i-get-the-behavior-of-gnus-readlink-f-on-a-mac
ADATD=$(python -c 'import os,sys;print(os.path.realpath(sys.argv[1]))' $DATD)
AOUTD=$(python -c 'import os,sys;print(os.path.realpath(sys.argv[1]))' $OUTD)

# /output in container maps to $OUTD on host
ARG="-o /output"
CMD="bash $PROCESS $@ $ARG $REF $BAM"
DCMD="docker run -v $ADATD:/data -v $AOUTD:/output -it $IMAGE $CMD"
>&2 echo Running: $DCMD
eval $DCMD


