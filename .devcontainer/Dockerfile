FROM debian:trixie-slim

# For intel-mkl
RUN sed -i 's/^Components: main.*$/Components: main contrib non-free non-free-firmware/g' /etc/apt/sources.list.d/debian.sources

# update packages and install dev dependencies (python3, pytest)
RUN apt-get update && apt-get install -y cmake python3 python3-pip python-is-python3

# Install mkl, openmpi
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y intel-mkl openmpi-bin libopenmpi-dev

# Install sudo for docker user
RUN apt-get install -y sudo
# Create docker user
ENV UNAME=docker GID=1000 UID=1000
RUN groupadd -g $GID -o $UNAME && useradd -m -u $UID -g $GID -G sudo -o -s /bin/bash $UNAME && \
    echo "$UNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Intel MKL
RUN apt-get install -y wget gpg && \
mkdir -p /usr/share/keyrings && \
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | \
gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null && \
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | \
tee /etc/apt/sources.list.d/oneAPI.list && \
apt-get update && \
apt-get install -y intel-oneapi-mkl

# Intel oneapi
# RUN apt-get install -y intel-oneapi-compiler-fortran intel-oneapi-openmp intel-oneapi-mpi intel-oneapi-mpi-devel

# Change to docker user
USER ${UNAME}
WORKDIR /workspace

# Add intel oneapi to path
RUN echo "source /opt/intel/oneapi/setvars.sh --force > /dev/null" >> ~/.bashrc
# Install dependencies
RUN pip install --break-system-packages black flake8 fortls fprettify pytest
