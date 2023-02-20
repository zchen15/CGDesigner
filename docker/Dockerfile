FROM ubuntu
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y
RUN apt-get install -y wget diffutils grep unzip python3-pip
RUN apt-get install -y vim git
RUN apt-get install -y openssh-server openssh-client
RUN apt-get install -y tmux
RUN pip3 install boto3
RUN pip3 install biopython
RUN pip3 install numpy
RUN pip3 install pandas
RUN pip3 install matplotlib
RUN pip3 install scikit-learn
RUN pip3 install scipy
RUN pip3 install psutil
RUN pip3 install kubernetes
RUN pip3 install ipython
COPY nupack-latest /root/nupack-latest
RUN pip3 install -U nupack -f /root/nupack-latest/package
COPY . /root/
WORKDIR /root
RUN pip3 install .
