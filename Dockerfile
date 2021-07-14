# syntax=docker/dockerfile:1

FROM ubuntu:20.04

WORKDIR /app

COPY requirements.txt requirements.txt

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt install -y python3.9
RUN apt install -y python3-pip
RUN pip3 install -r requirements.txt

COPY commands.sh /scripts/commands.sh

RUN ["chmod", "+x", "/scripts/commands.sh"]

COPY . .

CMD ["/scripts/commands.sh"]


