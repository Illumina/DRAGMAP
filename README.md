# Dragmap 

Dragmap is the Dragen mapper/aligner Open Source Software.

## Installation

### Prerequisites

Compilation was tested on CentOS 7

* C++11 compatible compiler (e.g gcc-c++ >= 4.8.5-36.el7_6.2)
* GNU make >= 3.82
* Boost library :  boost169-devel >= 1.69.0-1.el7
* For unit tests : googletest (>= v1.6)
* Hardware: x86_86, 64GB RAM minimum
* OS: Centos >= 7.7

### Install


The basic procedure is

    make

Binary will be generated in ./build/release/


Then optionally, to install to /usr/bin/

    make install



By default make will compile and launch unit tests. To disable unit tests, use HAS_GTEST=0, e.g. :


    HAS_GTEST=0 make


To compile with unit tests, if google test was installed in user space, it might be required to set GTEST_ROOT and LD_LIBRARY_PATH to where gtest was installed, e.g. : 

    export GTEST_ROOT=/home/username/lib/gtest
    export LD_LIBRARY_PATH=/home/username/lib/gtest/lib




#### Other variables controlling the build process:


* GCC_BASE 
* CXX 
* BOOST_ROOT 
* BOOST_INCLUDEDIR 
* BOOST_LIBRARYDIR 





## Basic command line usage 

### Command line options 

    dragen-os --help


### Build hash table of a reference fasta file 

    dragen-os --build-hash-table true --ht-reference reference.fasta  --output-directory /home/data/reference/


### Align paired-end reads :

Output result to standard output 

    dragen-os -r /home/data/reference/ -1 reads_1.fastq.gz -2 reads_2.fastq.gz >  result.sam

Or directly to a file :

    dragen-os -r /home/data/reference/ -1 reads_1.fastq.gz -2 reads_2.fastq.gz --output-directory /home/data/  --output-file-prefix result

### Align single-end reads :

    dragen-os -r /home/data/reference/ -1 reads_1.fastq.gz  >  result.sam

