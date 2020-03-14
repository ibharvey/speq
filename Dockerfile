FROM ubuntu:bionic AS build

RUN apt-get update && \
	apt-get install -y build-essential git cmake autoconf libtool pkg-config zlib1g-dev libbz2-dev

WORKDIR /speq

COPY include/ /speq/include/
COPY src/ /speq/src/
COPY deps/ /speq/deps/
COPY CMakeLists.txt /speq/

WORKDIR /speq/build

RUN cmake .. && make

FROM ubuntu:bionic AS release

WORKDIR /bin

COPY --from=build /speq/build/src/speq ./

ENTRYPOINT ["./speq"]
