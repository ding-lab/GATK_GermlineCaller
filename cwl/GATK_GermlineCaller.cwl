class: CommandLineTool
cwlVersion: v1.0
id: somatic_sv_workflow
baseCommand:
  - /bin/bash
  - /opt/GATK_GermlineCaller/src/process_sample_parallel.sh
inputs:
  - id: bam
    type: File
    inputBinding:
      position: 1
    label: Input BAM/CRAM
    secondaryFiles: ${if (self.nameext === ".bam") {return self.basename + ".bai"} else {return self.basename + ".crai"}}
  - id: reference
    type: File
    inputBinding:
      position: 2
    label: Reference FASTA
    secondaryFiles:
      - .fai
  - id: dryrun
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-d'
    label: dry run
    doc: 'Print out commands but do not execute, for testing only'
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
    dockerPull: 'mwyczalkowski/gatk_germlinecaller'
  - class: InlineJavascriptRequirement
