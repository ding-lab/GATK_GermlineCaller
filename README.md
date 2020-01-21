# GATK_GermlineCaller

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

Parameters

germline_variant_snakemake
gatk HaplotypeCaller -R genome_fa -I bam -L interval -O output --standard-min-confidence-threshold-for-calling 30.0

somaticswrapper
gatk HaplotypeCaller -I bam -O output -R genome_fa -L interval -RF NotDuplicateReadFilter -RF MappingQualityReadFilter -RF MappedReadFilter
