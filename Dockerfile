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

COPY . .

CMD ["python3", "Main.py", "fichier/COL1A1_trypsin.csv", "results/"]


