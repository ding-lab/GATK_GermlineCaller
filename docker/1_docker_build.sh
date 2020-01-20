IMAGE="mwyczalkowski/GATK_GermlineCaller"

cd ..
docker build -t $IMAGE -f docker/Dockerfile .
