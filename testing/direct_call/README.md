Testing of basic functionality running on test dataset
Before running, uncompress reference in ../demo_data with,
	tar -xvjf Homo_sapiens_assembly19.COST16011_region.fa.tar.bz2


Test as, (TODO confirm filenames)

bash 1_run_direct.sh -d
	Dry run to test basic configuration

bash 1_run_direct.sh 
	This runs GATK HaplotypeCaller on entire dataset, writes output/XXX and output/YYY

bash 2_run_direct_parallel.sh
	Runs GATK HaplotypeCaller on four regions in parallel, then merges them to write final output files
	output.parallel/XXX and output.parallel/YYY
