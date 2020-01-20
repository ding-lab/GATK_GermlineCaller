# Basic script to start docker container

source docker_image.sh

# Starting with no mounted volumes:
docker run -it $IMAGE
