FROM debian:bullseye-slim

# For intel-mkl
RUN echo "deb http://deb.debian.org/debian buster main contrib non-free \
    deb-src http://deb.debian.org/debian buster main contrib non-free \
    deb http://deb.debian.org/debian-security/ buster/updates main contrib non-free \
    deb-src http://deb.debian.org/debian-security/ buster/updates main contrib non-free \
    deb http://deb.debian.org/debian buster-updates main contrib non-free \
    deb-src http://deb.debian.org/debian buster-updates main contrib non-free" >> /etc/apt/sources.list

# update packages and install dev dependencies (python3, pytest)
RUN apt-get update && apt-get install -y cmake python3 python3-pip python-is-python3 && pip3 install pytest

# Install mkl, openmpi
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y intel-mkl openmpi-bin libopenmpi-dev

# Install sudo for docker user
RUN apt-get install -y sudo
# Create docker user
ENV UNAME=docker GID=1000 UID=1000
RUN groupadd -g $GID -o $UNAME && useradd -m -u $UID -g $GID -G sudo -o -s /bin/bash $UNAME && \
    echo "$UNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Intel MKL
RUN apt-get install -y wget && \
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list && \
apt-get update && \
apt-get install -y intel-oneapi-mkl

# Change to docker user
USER ${UNAME}
WORKDIR /workspace

RUN echo "source /opt/intel/oneapi/setvars.sh --force > /dev/null" >> ~/.bashrc

# Intel oneapi
# RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
#     apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
#     rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB && \
#     echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list && \
#     apt-get update && \
#     apt-get install -y intel-oneapi-compiler-fortran intel-oneapi-openmp intel-oneapi-mpi intel-oneapi-mpi-devel && \
#     echo "source /opt/intel/oneapi/setvars.sh" >> ~/.bashrc
