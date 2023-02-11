FROM ubuntu

WORKDIR /home/

COPY htslib-1.16.tar.bz2 .
COPY build_htslib.sh .

RUN apt-get update 
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN apt-get install -y unzip build-essential zlib1g-dev autoconf libbz2-dev liblzma-dev libcurl4-openssl-dev cmake

RUN ./build_htslib.sh

COPY CMakeLists.txt ./
COPY *.h ./
COPY *.cpp ./
ADD libs ./libs

RUN cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo . && make

RUN apt-get install -y python python-dev curl
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
RUN python get-pip.py

RUN pip install pysam pyfaidx numpy

COPY random_pos_generator.py surveyor.py ./

ENTRYPOINT [ "/usr/bin/python", "/home/surveyor.py" ]
