# GATK_GermlineCaller

Call germline variants using GATK4 HaplotypeCaller
Can operate on multiple regions with a passed CHRLIST file

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

## CWL 

IMAGE="mwyczalkowski/gatk_germlinecaller:20200608"

## CHRLIST

CHRLIST is a file which can take arbitrary genomic regions in a format accepted by GATK HaplotypeCaller.
Generally, a listing of all chromosomes will suffice

## Testing

`./testing` directory has demo data which can be quickly used to exercise different parts of pipeline
Pipeline can be called in 3 contexts:
* Direct, but entering docker container and running from command line 
* Docker, by invoking a docker run with the requested command
* CWL, using CWL workflow manager
  * Rabix and cromwell are supported

## Production

Setting `finalize` parameter to `true` will compress all intermediate files and logs
It also deletes intermediate raw vcf file
Setting `do_index` to true will compress output VCF files and create .csi index,
  and also remove the original (uncompressed) .vcf and .vcf.idx files
For production runs setting both to true is recommended

## Parameters

Not clear what additional parameters should be used for HaplotypeCaller.
germline_variant_snakemake has the argument 
    gatk HaplotypeCaller --standard-min-confidence-threshold-for-calling 30.0

## Author

Matthew Wyczalkowski <m.wyczalkowski@wustl.edu>

