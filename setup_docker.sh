
# Is the ubuntu:seqan3 image already in the docker daemon?
IID="$(docker images -q ubuntu:seqan3 2> /dev/null)"

# Build the ubuntu:seqan3 docker image
docker build -f docker/cmake.dkr -t ubuntu:seqan3 .

# If the ubuntu:seqan3 image was pre-built
# then there would now be a dangling image
# so we will remove that old image.
if [[ "${IID}" != "" ]]; then
	docker rmi "${IID}"	
fi

# Install the speq program itself
docker build -f docker/install.dkr -t speq .

