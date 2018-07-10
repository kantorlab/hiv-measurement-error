FROM ubuntu:16.04
LABEL maintainer "Mark Howison <mhowison@brown.edu>"
LABEL repository kantorlab
LABEL image hiv-measurement-error
LABEL tag v1

RUN apt-get update -y
RUN apt-get install -y bc bzip2 git openjdk-9-jre-headless unzip util-linux wget

RUN wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh \
 && bash Anaconda3-5.2.0-Linux-x86_64.sh -b \
 && rm Anaconda3-5.2.0-Linux-x86_64.sh

ENV PATH /root/anaconda3/bin:$PATH

RUN conda install -y -c kantorlab blastn=2.7.1 hivmmer=0.1.2 iva=1.0.9 mafft=7.313 matplotlib quasitools=0.3.1 scons=3.0.1.1 sra-tools=2.9.1.1 trimmomatic=0.36 seaborn
RUN conda clean -ay

RUN apt-get autoremove -y
RUN apt-get clean -y
