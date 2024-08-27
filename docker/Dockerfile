# Use the official Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Update the package list and install the required packages
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    git \
    gcc \
    gfortran \
    libopenmpi-dev \
    netcdf-bin \
    openmpi-bin \
    linux-tools-common \
    libnetcdf-dev \
    libnetcdff-dev \
    ca-certificates \
    autoconf \
    build-essential

# Clone and compile MOM6
RUN git clone https://github.com/NOAA-GFDL/MOM6-examples.git && \
    cd MOM6-examples && \
    git submodule update --init --recursive src/MOM6 && \
    cd src/MOM6/ac/deps && make -j && \
    cd .. && \
    autoreconf && \
    cd .. && \
    mkdir -p build && \
    cd build && \
    ../ac/configure && \
    make -j

# Get the MOM6 cases and setup the double gyre example
RUN mkdir -p MOM6-examples/ocean_only/double_gyre/RESTART

WORKDIR /MOM6-examples/ocean_only/double_gyre
CMD [ "mpirun", "--allow-run-as-root", "-n", "4", "/MOM6-examples/src/MOM6/build/MOM6" ]




