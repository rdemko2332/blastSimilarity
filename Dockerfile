FROM ubuntu:22.04

Label maintainer="rdemko2332@gmail.com"

WORKDIR /usr/bin/

RUN apt-get -qq update --fix-missing && \
    apt-get install -y \
    wget \
    perl \
    libgomp1 

# Setting up ncbi blast tools
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz \
  && tar -zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz \
  && rm -rf ncbi-blast-2.13.0+-x64-linux.tar.gz

# Adding Perl module files and fixZip.pl to usr/bin/
ADD /lib/perl/*.pm /lib/perl/
ADD /bin/*.pl /usr/bin/

# Making all blast tools executable
RUN chmod +x * \
  && cd ncbi-blast-2.13.0+/bin  \
  &&  chmod +x *

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/bin/ncbi-blast-2.13.0+/bin:/usr/bin/ncbi-blast-2.13.0+:/lib/perl/

ENV PERL5LIB=/lib/perl/

WORKDIR /work


