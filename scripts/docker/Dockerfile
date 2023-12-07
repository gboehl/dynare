# Copyright 2023 Dynare Team
# This file is part of Dynare.
#
# Dynare is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dynare is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

##############################################################################################################
# INFORMATION:                                                                                               #
# This Dockerfile installs Dynare (https://dynare.org)                                                       #
#  - with MATLAB including the following toolboxes (https://github.com/mathworks-ref-arch/matlab-dockerfile) #
#    MATLAB Symbolic_Math_Toolbox Statistics_and_Machine_Learning_Toolbox Optimization_Toolbox               #
#    Econometrics_Toolbox Parallel_Computing_Toolbox Control_System_Toolbox Global_Optimization_Toolbox      #
#  - with Octave including the following toolboxes (using debian packages):                                  #
#    octave-control octave-econometrics octave-io octave-statistics octave-struct octave-parallel            #
# It also renames the GCC libraries of MATLAB (similar to the matlab-support package) to avoid conflicts     #
# with the system libraries and adds the path of Dynare to the MATLAB and Octave startup scripts.            #
#                                                                                                            #
# MATLAB LICENSE:                                                                                            #
# The container is created using a network license, so no information on the license is inside the container #
# see https://git.dynare.org/dynare/dynare/docker/README.md#matlab-license for more information.             #
##############################################################################################################

# Default values which MATLAB and Dynare release to install in the container
# The Dynare release must conform to a corresponding tag on https://git.dynare.org/dynare/dynare
# MATLAB release must conform to a corresponding tag on https://hub.docker.com/r/mathworks/matlab/tags
# Octave version the one shipped with the Ubuntu version used in the base container (in 20.04 it is 6.4.0)
ARG MATLAB_RELEASE=R2023a
ARG DYNARE_RELEASE=5.4

# Specify the list of products to install into MATLAB with mpm
ARG MATLAB_PRODUCT_LIST="Symbolic_Math_Toolbox Statistics_and_Machine_Learning_Toolbox Optimization_Toolbox Econometrics_Toolbox Parallel_Computing_Toolbox Control_System_Toolbox Global_Optimization_Toolbox"

# Specify MATLAB Install Location.
ARG MATLAB_INSTALL_LOCATION="/opt/matlab/${MATLAB_RELEASE}"

# Specify license server information using the format: port@hostname 
ARG LICENSE_SERVER

# Specify the base image with MATLAB installed.
FROM mathworks/matlab:${MATLAB_RELEASE}
USER root

# Declare build arguments to use at the current build stage.
ARG MATLAB_RELEASE
ARG MATLAB_PRODUCT_LIST
ARG MATLAB_INSTALL_LOCATION
ARG LICENSE_SERVER
ARG DYNARE_RELEASE

# Install mpm dependencies.
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get install --no-install-recommends --yes \
    wget \
    unzip \
    ca-certificates \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Run mpm to install additional toolboxes for MATLAB in the target location and delete the mpm installation afterwards.
# If mpm fails to install successfully, then print the logfile in the terminal, otherwise clean up.
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm \ 
    && chmod +x mpm \
    && ./mpm install \
    --release=${MATLAB_RELEASE} \
    --destination=${MATLAB_INSTALL_LOCATION} \
    --products ${MATLAB_PRODUCT_LIST} \
    || (echo "MPM Installation Failure. See below for more information:" && cat /tmp/mathworks_root.log && false) \
    && rm -f mpm /tmp/mathworks_root.log

# Install dynare dependencies.
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get install --no-install-recommends --yes \
    build-essential \
    gfortran \
    libboost-graph-dev \
    libgsl-dev \
    libmatio-dev \
    libslicot-dev \
    libslicot-pic \
    libsuitesparse-dev \
    flex \
    libfl-dev \
    bison \
    autoconf \
    automake \
    texlive \
    texlive-publishers \
    texlive-latex-extra \
    texlive-fonts-extra \
    texlive-font-utils \
    texlive-latex-recommended \
    texlive-science \
    texlive-plain-generic \
    lmodern \
    python3-sphinx \
    tex-gyre \
    latexmk \
    libjs-mathjax \
    doxygen \
    x13as \
    liboctave-dev \
    octave-control \
    octave-econometrics \
    octave-io \
    octave-statistics \
    octave-struct \
    octave-parallel \
    gnuplot \
    fonts-freefont-otf \
    ghostscript \
    epstool \
    git \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Rename libraries (see matlab-support package: https://salsa.debian.org/debian/matlab-support/-/blob/master/debian/matlab-support.postinst)
RUN if [ -f "${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libgcc_s.so.1" ]; then \
        mv ${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libgcc_s.so.1 ${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libgcc_s.so.1.bak; \
    fi && \
    if [ -f "${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libstdc++.so.6" ]; then \
        mv ${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libstdc++.so.6 ${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libstdc++.so.6.bak; \
    fi && \
    if [ -f "${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libgfortran.so.5" ]; then \
        mv ${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libgfortran.so.5 ${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libgfortran.so.5.bak; \
    fi && \
    if [ -f "${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libquadmath.so.0" ]; then \
        mv ${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libquadmath.so.0 ${MATLAB_INSTALL_LOCATION}/sys/os/glnxa64/libquadmath.so.0.bak; \
    fi && \
    if [ -f "${MATLAB_INSTALL_LOCATION}/bin/glnxa64/libfreetype.so.6" ]; then \
        mv ${MATLAB_INSTALL_LOCATION}/bin/glnxa64/libfreetype.so.6 ${MATLAB_INSTALL_LOCATION}/bin/glnxa64/libfreetype.so.6.bak; \
    fi
# Fix for epstopdf latex errors, i.e. LIBTIFF_4.0 not found
RUN if [ -f "${MATLAB_INSTALL_LOCATION}/bin/glnxa64/libtiff.so.5" ]; then \
        mv ${MATLAB_INSTALL_LOCATION}/bin/glnxa64/libtiff.so.5 ${MATLAB_INSTALL_LOCATION}/bin/glnxa64/libtiff.so.5.bak ; \
    fi

# Note: Uncomment one of the following two ways to configure the license server.
# WE DO NOT WANT OPTION 2!!!!
# Option 1. Specify the host and port of the machine that serves the network licenses
# if you want to store the license information in an environment variable. This
# is the preferred option. You can either use a build variable, like this: 
# --build-arg LICENSE_SERVER=27000@MyServerName or you can specify the license server 
# directly using: ENV MLM_LICENSE_FILE=27000@flexlm-server-name
ENV MLM_LICENSE_FILE=$LICENSE_SERVER
# Option 2. Alternatively, you can put a license file into the container (WE DON'T WANT THIS!!!)
# Enter the details of the license server in this file and uncomment the following line.
# COPY network.lic ${MATLAB_INSTALL_LOCATION}/licenses/

# Get Dynare sources as matlab user
USER matlab
WORKDIR /home/matlab
RUN git clone --depth 1 --branch ${DYNARE_RELEASE} --recurse-submodules https://git.dynare.org/dynare/dynare.git

# Compile Dynare
USER matlab
WORKDIR /home/matlab
RUN cd dynare \
    && autoreconf -si \
    && ./configure --with-matlab=${MATLAB_INSTALL_LOCATION} MATLAB_VERSION=${MATLAB_RELEASE} \
    && make -j$(($(nproc)+1))

# Add path of dynare to startup script for Octave.
RUN echo "addpath /home/matlab/dynare/matlab" >> /home/matlab/.octaverc

# Add path of dynare to startup script for MATLAB.
# Note that if startup.m file exists (in newer MATLAB containers), it is a MATLAB function
# and the last line is an "end", so we append the path to the second-to-last line
# For some reason we have to do this as root, otherwise the file is not writable
USER root
RUN filename="/home/matlab/Documents/MATLAB/startup.m" && \
    if [ -f ${filename} ]; then \
      tempfile=$(mktemp) && \
      total_lines=$(wc -l < "$filename") && \
      line_number=$((total_lines - 1)) && \
      head -n "$line_number" "$filename" > "$tempfile" && \
      echo "addpath('/home/matlab/dynare/matlab');" >> "$tempfile" && \
      tail -n 1 "$filename" >> "$tempfile" && \
      mv "$tempfile" "$filename" ; \
    else \
      echo "addpath('/home/matlab/dynare/matlab');" >> "$filename" ; \
    fi && \
    chown matlab:matlab "$filename"

# Set user and work directory.
USER matlab
WORKDIR /home/matlab