class: CommandLineTool
cwlVersion: v1.0
id: GATK_GermlineCaller
baseCommand:
  - /bin/bash
  - /opt/GATK_GermlineCaller/src/process_sample_parallel.sh
inputs:
  - id: reference
    type: File
    inputBinding:
      position: 1
    label: Reference FASTA
    secondaryFiles:
      - .fai
      - ^.dict
  - id: bam     # CRAM support untested
    type: File
    inputBinding:
      position: 2
    label: Input BAM
    secondaryFiles: ${if (self.nameext === ".bam") {return self.basename + ".bai"} else {return self.basename + ".crai"}}
  - id: chrlist
    type: File?
    inputBinding:
      position: 0
      prefix: '-c'
    doc: List of genomic regions
    label: Genomic regions
  - id: njobs
    type: int?
    inputBinding:
      position: 0
      prefix: '-j'
    label: N parallel jobs
    doc: 'Number of jobs to run in parallel mode'
  - id: dryrun
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-d'
    label: dry run
    doc: 'Print out commands but do not execute, for testing only'
  - id: finalize
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-F'
    label: finalize
    doc: 'Compress intermediate data and logs, remove raw VCF'
  - id: compress_output
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-I'
    label: Compress output
    doc: 'Compress and index output VCF files'
  - id: HC_ARGS
    type: string?
    inputBinding:
      position: 0
      prefix: '-C'
    label: HaplotypeCaller args
    doc: GATK HaplotypeCaller arguments
  - id: SV_SNP_ARGS
    type: string?
    inputBinding:
      position: 0
      prefix: '-R'
    label: SelectVariants SNP args
    doc: SelectVariants SNP arguments
  - id: SV_INDEL_ARGS
    type: string?
    inputBinding:
      position: 0
      prefix: '-S'
    doc: SelectVariants INDEL arguments
    label: SelectVariants INDEL args
outputs:
  - id: snp_vcf
    type: File?
    outputBinding:
      glob: ${if (inputs.compress_output ) {return "output/GATK.snp.Final.vcf.gz" } else {return "output/GATK.snp.Final.vcf"}}
  - id: indel_vcf
    type: File?
    outputBinding:
      glob: ${if (inputs.compress_output ) {return "output/GATK.indel.Final.vcf.gz" } else {return "output/GATK.indel.Final.vcf"}}
label: GATK_GermlineCaller
requirements:
  - class: ResourceRequirement
    ramMin: 12000
  - class: DockerRequirement
    dockerPull: 'mwyczalkowski/gatk_germlinecaller:20200608'
  - class: InlineJavascriptRequirement
