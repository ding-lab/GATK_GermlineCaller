Testing on compute1

Based on part on TinDaisy testing on compute1 here: /home/m.wyczalkowski/Projects/TinDaisy/TinDaisy/demo/MutectDemo/compute1-dev

Currently (1/22/20), test run fails with WORKFLOW_ROOT="/data/Active/cromwell-data"

Changing WORKFLOW_ROOT="/home/m.wyczalkowski/Projects/dat/cromwell-data"
-> Succeeds

Successful run of 30
/home/m.wyczalkowski/Projects/dat/cromwell-data/cromwell-workdir/cromwell-executions/GATK_GermlineCaller.cwl/d1534412-a8b4-4c01-87d5-a4704aa51442/call-GATK_GermlineCaller.cwl/execution/output/GATK.indel.Final.vcf

Successful run of 35
/home/m.wyczalkowski/Projects/dat/cromwell-data/cromwell-workdir/cromwell-executions/GATK_GermlineCaller.cwl/32814174-62c6-410c-81dc-f73898f0ec59/call-GATK_GermlineCaller.cwl/execution/output/GATK.indel.Final.vcf

WID="d1534412-a8b4-4c01-87d5-a4704aa51442"
curl -k -s -X GET http://localhost:8000/api/workflows/v1/$WID/status -H "accept: application/json"
{"status":"Succeeded","id":"32814174-62c6-410c-81dc-f73898f0ec59"}
-> cq will work!
