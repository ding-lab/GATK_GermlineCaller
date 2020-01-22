class: CommandLineTool
cwlVersion: v1.0
id: somatic_sv_workflow
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
  - id: bam
    type: File
    inputBinding:
      position: 2
    label: Input BAM/CRAM
    secondaryFiles: ${if (self.nameext === ".bam") {return self.basename + ".bai"} else {return self.basename + ".crai"}}
  - id: chrlist
    type: File?
    inputBinding:
      position: 0
      prefix: '-c'
    label: List of genomic regions
  - id: njobs
    type: int?
    inputBinding:
      position: 0
      prefix: '-j'
    label: Parallel job count
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
    doc: 'Compress intermediate data and logs'
  - id: HC_ARGS
    type: string?
    inputBinding:
      position: 0
      prefix: '-C'
    label: GATK HaplotypeCaller arguments
  - id: SV_SNP_ARGS
    type: string?
    inputBinding:
      position: 0
      prefix: '-R'
    label: SelectVariants SNP arguments
  - id: SV_INDEL_ARGS
    type: string?
    inputBinding:
      position: 0
      prefix: '-S'
    label: SelectVariants INDEL arguments
outputs:
  - id: snp_vcf
    type: File?
    outputBinding:
      glob: output/GATK.snp.Final.vcf
  - id: indel_vcf
    type: File?
    outputBinding:
      glob: output/GATK.indel.Final.vcf
label: GATK_GermlineCaller
requirements:
  - class: ResourceRequirement
    ramMin: 8000
  - class: DockerRequirement
    dockerPull: mwyczalkowski/gatk_germlinecaller
  - class: InlineJavascriptRequirement
