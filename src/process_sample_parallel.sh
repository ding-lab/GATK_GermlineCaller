#!/bin/bash

# Matthew Wyczalkowski <m.wyczalkowski@wustl.edu>
# https://dinglab.wustl.edu/

read -r -d '' USAGE <<'EOF'
Run GATK germline caller, possibly for multiple intervals in parallel,
and generate one VCF file for SNVs and one for INDELs

Usage: 
  process_sample_parallel.sh [options] REF BAM

Output: 
    OUTD/GATK.snp.Final.vcf
    OUTD/GATK.indel.Final.vcf

Options:
-h : print usage information
-d : dry-run. Print commands but do not execute them
-1 : stop after iteration over CHRLIST
-c CHRLIST: File listing genomic intervals over which to operate
-j JOBS: if parallel run, number of jobs to run at any one time.  If 0, run sequentially.  Default: 4
-o OUTD: set output root directory.  Default ./output
-F : finalize run by compressing per-region output and logs

The following arguments are passed to process_sample.sh directly:
-C HC_ARGS : pass args to `gatk HaplotypeCaller`
-R SV_SNP_ARGS : pass SV_SNP_ARGS to `gatk SelectVariants -select-type SNP -select-type MNP`
-S SV_INDEL_ARGS : pass SV_INDEL_ARGS to `gatk SelectVariants -select-type INDEL`

For single region, calls look like,:
  gatk HaplotypeCaller -R REF -I BAM 
  gatk SelectVariants -O GATK.snp.Final.vcf -select-type SNP -select-type MNP 
  gatk SelectVariants -O GATK.indel.Final.vcf -select-type INDEL

For multiple regions (specified by -c CHRLIST), calls are like,
  for CHR in CHRLIST
    gatk HaplotypeCaller -R REF -I BAM -L CHR
    gatk SelectVariants -O CHR_SNP -select-type SNP -select-type MNP 
    gatk SelectVariants -O CHR_INDEL -select-type INDEL
  bcftools concat -o GATK.snp.Final.vcf
  bcftools concat -o GATK.indel.Final.vcf

CHRLIST is a file listing genomic intervals over which to operate, with each
line passed to `gatk HaplotypeCaller -L`. 

In general, if CHRLIST is defined, jobs will be submitted in parallel mode: use
GNU parallel to loop across all entries in CHRLIST, running -j JOBS at a time,
and wait until all jobs completed.  Output logs written to OUTD/logs/GATK.$CHR.log
Parallel mode can be disabled with -j 0.

EOF

source /opt/GATK_GermlineCaller/src/utils.sh
SCRIPT=$(basename $0)

# Background on `parallel` and details about blocking / semaphores here:
#    O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
#    ;login: The USENIX Magazine, February 2011:42-47.
# [ https://www.usenix.org/system/files/login/articles/105438-Tange.pdf ]

# set defaults
NJOBS=4
DO_PARALLEL=0
OUTD="./output"
PROCESS="/opt/GATK_GermlineCaller/src/process_sample.sh"
BCFTOOLS="/opt/miniconda/bin/bcftools"

# http://wiki.bash-hackers.org/howto/getopts_tutorial
while getopts ":hd1c:j:o:C:R:S:F" opt; do
  case $opt in
    h)
      echo "$USAGE"
      exit 0
      ;;
    d)  # example of binary argument
      >&2 echo "Dry run" 
      DRYRUN=1
      ;;
    1) 
      JUSTONE=1
      ;;
    c) 
      CHRLIST_FN=$OPTARG
      DO_PARALLEL=1
      ;;
    j) 
      NJOBS=$OPTARG  
      ;;
    o) 
      OUTD=$OPTARG
      ;;
    C) 
      PS_ARGS="$PS_ARGS -C \"$OPTARG\""
      ;;
    R) 
      PS_ARGS="$PS_ARGS -R \"$OPTARG\""
      ;;
    S) 
      PS_ARGS="$PS_ARGS -S \"$OPTARG\""
      ;;
    F) 
      FINALIZE=1
      ;;
    \?)
      >&2 echo "$SCRIPT: ERROR: Invalid option: -$OPTARG"
      >&2 echo "$USAGE"
      exit 1
      ;;
    :)
      >&2 echo "$SCRIPT: ERROR: Option -$OPTARG requires an argument."
      >&2 echo "$USAGE"
      exit 1
      ;;
  esac
done
shift $((OPTIND-1))

if [ "$#" -ne 2 ]; then
    >&2 echo ERROR: Wrong number of arguments
    >&2 echo "$USAGE"
    exit 1
fi

REF=$1;   confirm $REF
BAM=$2;   confirm $BAM

# Read CHRLIST_FN to get list of elements
# These were traditionally individual chromosomes, 
# but can be other regions as well
if [ $CHRLIST_FN ]; then
    confirm $CHRLIST_FN
    CHRLIST=$(cat $CHRLIST_FN)
    # Will need to merge multiple VCFs
else
    # Will not need to merge multiple VCFs
    CHRLIST="Final"
    NO_CHRLIST=1
fi
    
# Output, tmp, and log files go here
mkdir -p $OUTD

# Per-region output goes here
if [ ! "$NO_CHRLIST" ]; then
    OUTDR="$OUTD/regions"
    mkdir -p $OUTDR
fi

# CHRLIST newline-separated list of regions passed to GATK HaplotypeCaller -L 
LOGD="$OUTD/logs"
mkdir -p $LOGD

NOW=$(date)
MYID=$(date +%Y%m%d%H%M%S)

if [ $NJOBS == "0" ]; then 
    DO_PARALLEL=0
fi

if [ $DO_PARALLEL == 1 ]; then
    >&2 echo [ $NOW ]: Parallel run 
    >&2 echo . 	  Looping over $CHRLIST
    >&2 echo . 	  Parallel jobs: $NJOBS
    >&2 echo . 	  Log files: $LOGD
else
    >&2 echo [ $NOW ]: Single region at a time
    >&2 echo . 	  Looping over $CHRLIST
    >&2 echo . 	  Log files: $LOGD
fi

for CHR in $CHRLIST; do
    NOW=$(date)
    >&2 echo \[ $NOW \] : Processing $CHR

    STDOUT_FN="$LOGD/GATK_GermlineCaller.$CHR.out"
    STDERR_FN="$LOGD/GATK_GermlineCaller.$CHR.err"

    # core call to process_sample.sh
    if [ "$NO_CHRLIST" ]; then
        CMD="$PROCESS $PS_ARGS -o $OUTD -l Final $REF $BAM > $STDOUT_FN 2> $STDERR_FN"
    else
        CMD="$PROCESS $PS_ARGS -o $OUTDR -L $CHR $REF $BAM > $STDOUT_FN 2> $STDERR_FN"
    fi

    if [ $DO_PARALLEL == 1 ]; then
        JOBLOG="$LOGD/GATK_GermlineCaller.$CHR.log"
        CMD=$(echo "$CMD" | sed 's/"/\\"/g' )   # This will escape the quotes in $CMD
        CMD="parallel --semaphore -j$NJOBS --id $MYID --joblog $JOBLOG --tmpdir $LOGD \"$CMD\" "
    fi

    run_cmd "$CMD" $DRYRUN

    if [ "$JUSTONE" ]; then
        >&2 echo Exiting after one
        break
    fi
done

if [ $DO_PARALLEL == 1 ]; then
    # this will wait until all jobs completed
    CMD="parallel --semaphore --wait --id $MYID"
    run_cmd "$CMD" $DRYRUN
fi

# Now merge if we are looping over regions in CHRLIST
# Merged output will have same filename as if CHR were "Final"
# testing for globs from https://stackoverflow.com/questions/2937407/test-whether-a-glob-has-any-matches-in-bash
if [ ! "$NO_CHRLIST" ]; then

    # First merge the snp
    PATTERN="$OUTDR/GATK.snp.*.vcf"
    if stat -t $PATTERN >/dev/null 2>&1; then
        IN=`ls $PATTERN`
        OUT="$OUTD/GATK.snp.Final.vcf"
        CMD="$BCFTOOLS concat -o $OUT $IN"
        run_cmd "$CMD" $DRYRUN
        >&2 echo Final SNP output : $OUT
    else
        >&2 echo $SCRIPT : snp merge: no output found matching $PATTERN
    fi

    # then merge the indel
    PATTERN="$OUTDR/GATK.indel.*.vcf"
    if stat -t $PATTERN >/dev/null 2>&1; then
        IN=`ls $PATTERN`
        OUT="$OUTD/GATK.indel.Final.vcf"
        CMD="$BCFTOOLS concat -o $OUT $IN"
        run_cmd "$CMD" $DRYRUN
        >&2 echo Final INDEL output : $OUT
    else
        >&2 echo $SCRIPT : indel merge: no output found matching $PATTERN
    fi
fi

if [[ "$FINALIZE" ]] ; then

    LOGD="$OUTD/logs"
    TAR="$OUTD/logs.tar.gz"
    if [ -e $TAR ]; then
        >&2 echo WARNING: $TAR exists
        >&2 echo Skipping log finalize
    else
        CMD="tar -zcf $TAR $LOGD && rm -rf $LOGD"
        run_cmd "$CMD" $DRYRUN
        >&2 echo Logs in $LOGD is compressed as $TAR and deleted
    fi

    if [[ ! "$NO_CHRLIST" ]]; then
        TAR="$OUTD/regions.tar.gz"
        if [ -e $TAR ]; then
            >&2 echo WARNING: $TAR exists
            >&2 echo Skipping regions finalize
        else
            CMD="tar -zcf $TAR $OUTDR && rm -rf $OUTDR"
            run_cmd "$CMD" $DRYRUN
            >&2 echo Intermediate output in $OUTDR is compressed as $TAR and deleted
        fi
    fi
fi


NOW=$(date)
>&2 echo [ $NOW ] $SCRIPT : SUCCESS
