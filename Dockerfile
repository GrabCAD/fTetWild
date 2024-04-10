# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Avoid prompts from apt
ARG DEBIAN_FRONTEND=noninteractive

# Update and install necessary packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    wget \
    ca-certificates \
    jq \
    git \
    silversearcher-ag \
    ninja-build

RUN \
    rm -rf /var/lib/apt/lists/*

ADD . /usr/src/app

WORKDIR /usr/src/app

# fetch depends and setup build files
RUN cmake -B build_in_docker -S . -D CMAKE_BUILD_TYPE=Release -G Ninja

# build the dependencies and ftetwild
RUN cmake --build build_in_docker
