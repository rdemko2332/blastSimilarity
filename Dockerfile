FROM ubuntu:22.04
MAINTAINER rdemko2332@gmail.com
WORKDIR /work
RUN apt-get -qq update --fix-missing
RUN apt-get install -y wget perl libgomp1
RUN cd /usr/bin &&  wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz && tar -zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz && rm -rf ncbi-blast-2.13.0+-x64-linux.tar.gz
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/bin/ncbi-blast-2.13.0+/bin:/usr/bin/ncbi-blast-2.13.0+
RUN cd /usr/bin/ncbi-blast-2.13.0+/bin  &&  chmod +x *
COPY /bin/* /usr/bin/
RUN cd /usr/bin/ && chmod +x blastSimilarity
ENV PERL5LIB=/usr/bin/
COPY /data/* /work/
RUN makeblastdb -in db.fa -dbtype prot



