FROM ubuntu@sha256:51523b5adbc67853e73d7e5faff234614942f9ff3872f329d2bb59478baf53db
LABEL description="Builder for 3B on an ancient system"

# Since 10.04 (lucid) is long out of support, the packages have moved
RUN echo 'deb http://old-releases.ubuntu.com/ubuntu/ lucid main restricted universe' > /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y --force-yes openjdk-6-jre-headless && \
	apt-get install -y --force-yes openjdk-6-jdk wget zip vim && \
	apt-get install -y --force-yes libjpeg-dev libpng-dev libtiff-dev && \
	apt-get install -y --force-yes build-essential g++ 

RUN mkdir -p /tmp/deps /usr/local/lib /usr/local/include

WORKDIR /tmp/deps

#ADD ij145.zip   /tmp/deps
#ADD jpegsrc.v9a.tar.gz   /tmp/deps
#ADD libpng-1.6.18.tar.xz   /tmp/deps
#ADD tiff-3.9.7.tar.gz   /tmp/deps
#ADD zlib-1.2.8.tar.gz /tmp/deps

#Build lapack
ADD clapack.tgz   /tmp/deps
ADD clapack-make-linux.patch /tmp/deps
ADD clapack_mingw.patch /tmp/deps
WORKDIR /tmp/deps/CLAPACK-3.2.1
RUN cp make.inc.example make.inc && patch make.inc < ../clapack-make-linux.patch 
RUN patch -p1 < ../clapack_mingw.patch
RUN make -j8 blaslib && make -j8 f2clib && cd INSTALL && make ilaver.o slamch.o dlamch.o lsame.o && echo > second.c && cd .. && make -j8 lapacklib
RUN cp blas_LINUX.a /usr/local/lib/libblas.a && cp lapack_LINUX.a /usr/local/lib/liblapack.a && cp F2CLIBS/libf2c.a /usr/local/lib/libf2c.a

ADD TooN-2.0.tar.gz   /tmp/deps
WORKDIR /tmp/deps/TooN-2.0
RUN ./configure && make install


ADD gvars-3.0.tar.gz   /tmp/deps
WORKDIR /tmp/deps/gvars-3.0
RUN ./configure --without-head --without-lang && make -j8 && make install 

ADD libcvd-20121025.tar.gz   /tmp/deps
WORKDIR /tmp/deps/libcvd-20121025
RUN ./configure --disable-fast7 --disable-fast8 --disable-fast9 --disable-fast10 --disable-fast11 --disable-fast12 && make -j8 && make install

RUN mkdir -p /home/build
WORKDIR /home/build
