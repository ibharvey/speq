# Is the ubuntu:seqan3 image already in the docker daemon?
IID_SEQAN3="$(docker images -q ubuntu:seqan3 2> /dev/null)"

# Build the ubuntu:seqan3 docker image
docker build --iidfile iid.tmp -f ${1}/dockerfiles/ubuntu/cmake.dkr -t ubuntu:seqan3 ${1}

IID_SEQAN3_POST="$(cat iid.tmp)"
rm iid.tmp

# If the ubuntu:seqan3 image was pre-built
# then there would now be a dangling image
# so we will remove that old image.
if [[ "${IID_SEQAN3}" != "" ]]; then
        if [[ "${IID_SEQAN3}" != "${IID_SEQAN3_POST:7:12}" ]]; then
                docker rmi "${IID_SEQAN3}"
        fi
fi

# Is the speq:latest image already in the docker daemon?
IID_speq="$(docker images -q speq:latest 2> /dev/null)"

# Install the speq program itself
docker build --iidfile iid.tmp -f ${1}/dockerfiles/ubuntu/install.dkr -t speq ${1}

IID_speq_POST="$(cat iid.tmp)"
rm iid.tmp

if [[ "${IID_speq}" != "" ]]; then
        if [[ "${IID_speq}" != "${IID_speq_POST:7:12}" ]]; then
                docker rmi "${IID_speq}"
        fi
fi
