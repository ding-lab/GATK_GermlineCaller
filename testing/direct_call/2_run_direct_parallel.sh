# Example command to run within docker.  Typically, start docker first with 0_start_docker.sh

BAM="/opt/GATK_GermlineCaller/testing/demo_data/HCC1954.NORMAL.30x.compare.COST16011_region.bam"
REF="/opt/GATK_GermlineCaller/testing/demo_data/Homo_sapiens_assembly19.COST16011_region.fa"

PROCESS="/opt/GATK_GermlineCaller/src/process_sample_parallel.sh"

OUTD="output.parallel"
CHRLIST="chrlist.dat"

bash $PROCESS "$@" -c $CHRLIST  -o $OUTD $REF $BAM 

