source docker_image.sh
bsub -Is -q research-hpc -a "docker($IMAGE)" /bin/bash

