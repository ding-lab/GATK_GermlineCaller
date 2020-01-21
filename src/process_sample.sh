#/bin/bash

read -r -d '' USAGE <<'EOF'
Run GATK HaplotypeCaller and SelectVariants

Usage: process_sample.sh [options] reference.fa input.bam 
 
Options:
-h: Print this help message
-d : Dry run - output commands but do not execute them
-C HC_ARGS : pass args to `gatk HaplotypeCaller`
-R SV_SNP_ARGS : pass SV_SNP_ARGS to `gatk SelectVariants -select-type SNP -select-type MNP`
-S SV_INDEL_ARGS : pass SV_INDEL_ARGS to `gatk SelectVariants -select-type INDEL`
-L INPUT_INTERVAL : One or more genomic intervals over which to operate
  This is passed verbatim to HaplotypeCaller -L INPUT_INTERVAL
-l INTERVAL_LABEL : A short label for interval, used for filenames.  Default is INPUT_INTERVAL
-o OUTD : Output directory [ ./output ]

Output filenames:
    OUTD/GATK.snp.XXX.vcf
    OUTD/GATK.indel.XXX.vcf
where XXX is given by INTERVAL_LABEL
EOF

# Commands - templates:
# from germline_variant_snakemake:
#   gatk HaplotypeCaller -R {input.genome_fa} -I {input.bam} -L {input.interval} -O {output} --standard-min-confidence-threshold-for-calling 30.0
#   gatk SelectVariants -R {input.genome_fa} -V {input.input_vcf} -O {output} -select-type SNP -select-type MNP 
#   gatk SelectVariants -R {input.genome_fa} -V {input.input_vcf} -O {output} -select-type INDEL
#   bcftools concat -o   -- this does not take place here
# 
# From germlinewrapper:
#   gatk HaplotypeCaller  -I ${NBAM} -O ${rawvcf} -R $h38_REF -L $chr1 -RF NotDuplicateReadFilter -RF MappingQualityReadFilter -RF MappedReadFilter
#   gatk SelectVariants -R $h38_REF -V ${gvipvcf} -O ${snvvcf} -select-type SNP -select-type MNP
#   gatk SelectVariants -R $h38_REF -V ${gvipvcf} -O ${indelvcf} -select-type INDEL

# HaplotypeCaller documentation: https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
# SelectVariants documentation: https://gatk.broadinstitute.org/hc/en-us/articles/360037225432-SelectVariants

source /opt/GATK_GermlineCaller/src/utils.sh
SCRIPT=$(basename $0)

GATK="/gatk/gatk"

# Set defaults
OUTD="./output"
OUTVCF="final.SV.WGS.vcf"

# http://wiki.bash-hackers.org/howto/getopts_tutorial
while getopts ":hdC:R:S:L:l:o:" opt; do
  case $opt in
    h)
      echo "$USAGE"
      exit 0
      ;;
    d)  # binary argument
      DRYRUN=1
      ;;
    C) # value argument
      HC_ARGS="$HC_ARGS $OPTARG"
      ;;
    R) # value argument
      SV_SNP_ARGS="$SV_SNP_ARGS $OPTARG"
      ;;
    S) # value argument
      SV_INDEL_ARGS="$SV_INDEL_ARGS $OPTARG"
      ;;
    L) # value argument
      INPUT_INTERVAL="$OPTARG"
      ;;
    l) # value argument
      INTERVAL_LABEL="$OPTARG"
      ;;
    o) # value argument
      OUTD=$OPTARG
      ;;
    \?)
      >&2 echo "Invalid option: -$OPTARG"
      >&2 echo "$USAGE"
      exit 1
      ;;
    :)
      >&2 echo "Option -$OPTARG requires an argument."
      >&2 echo "$USAGE"
      exit 1
      ;;
  esac
done
shift $((OPTIND-1))


if [ "$#" -ne 2 ]; then
    >&2 echo Error: Wrong number of arguments
    >&2 echo "$USAGE"
    exit 1
fi

REF=$1
BAM=$2

confirm $BAM
confirm $REF

# IX forms part of suffix of output filename
# If INTERVAL_LABEL is given, IX takes that value
# otherwise, get value from INPUT_INTERVAL
if [ "$INTERVAL_LABEL" ]; then
    IX="$INTERVAL_LABEL"
else
    IX="$INPUT_INTERVAL"
fi

if [ "$INPUT_INTERVAL" ]; then
    HC_ARGS="$HC_ARGS -L $INPUT_INTERVAL"
fi

mkdir -p $OUTD
test_exit_status

# Run HaplotypeCaller
OUT1="$OUTD/GATK.raw.${IX}.vcf"
CMD1="$GATK HaplotypeCaller -R $REF -I $BAM -O $OUT1 $HC_ARGS"
run_cmd "$CMD1" $DRYRUN

# Run SelectVariants
OUT2A="$OUTD/GATK.snp.${IX}.vcf"
OUT2B="$OUTD/GATK.indel.${IX}.vcf"
CMD2A="$GATK SelectVariants -R $REF -V $OUT1 -O $OUT2A -select-type SNP -select-type MNP $SV_SNP_ARGS"
CMD2B="$GATK SelectVariants -R $REF -V $OUT1 -O $OUT2A -select-type INDEL $SV_INDEL_ARGS"
run_cmd "$CMD2A" $DRYRUN
run_cmd "$CMD2B" $DRYRUN

>&2 echo Written to $CMD2A and $CMD2B
>&2 echo $SCRIPT success.
