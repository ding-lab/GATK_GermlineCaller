Test data from manta distribution: https://github.com/Illumina/manta/tree/master/src/demo/data

Please uncompress reference with,
	tar -xvjf Homo_sapiens_assembly19.COST16011_region.fa.tar.bz2

java -jar /usr/picard/picard.jar CreateSequenceDictionary R=Homo_sapiens_assembly19.COST16011_region.fa O=Homo_sapiens_assembly19.COST16011_region.dict
