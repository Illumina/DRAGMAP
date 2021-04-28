# Dragen mapper/aligner - software version 

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



## Requirements

The binary was built using Centos 7. 
It requires the following dynamic libraries to run :

    ldd dragen-os 
	linux-vdso.so.1 =>  (0x00007ffe643f3000)
	libboost_system.so.1.53.0 => /lib64/libboost_system.so.1.53.0 (0x00007f9a7c88e000)
	libboost_filesystem.so.1.53.0 => /lib64/libboost_filesystem.so.1.53.0 (0x00007f9a7c677000)
	libboost_date_time.so.1.53.0 => /lib64/libboost_date_time.so.1.53.0 (0x00007f9a7c466000)
	libboost_thread-mt.so.1.53.0 => /lib64/libboost_thread-mt.so.1.53.0 (0x00007f9a7c24f000)
	libboost_system-mt.so.1.53.0 => /lib64/libboost_system-mt.so.1.53.0 (0x00007f9a7c04b000)
	libboost_iostreams.so.1.53.0 => /lib64/libboost_iostreams.so.1.53.0 (0x00007f9a7be31000)
	libboost_regex.so.1.53.0 => /lib64/libboost_regex.so.1.53.0 (0x00007f9a7bb2e000)
	libboost_program_options.so.1.53.0 => /lib64/libboost_program_options.so.1.53.0 (0x00007f9a7b8bc000)
	libz.so.1 => /lib64/libz.so.1 (0x00007f9a7b6a6000)
	libstdc++.so.6 => /lib64/libstdc++.so.6 (0x00007f9a7b39f000)
	librt.so.1 => /lib64/librt.so.1 (0x00007f9a7b197000)
	libgomp.so.1 => /lib64/libgomp.so.1 (0x00007f9a7af71000)
	libpthread.so.0 => /lib64/libpthread.so.0 (0x00007f9a7ad55000)
	libm.so.6 => /lib64/libm.so.6 (0x00007f9a7aa53000)
	libgcc_s.so.1 => /lib64/libgcc_s.so.1 (0x00007f9a7a83d000)
	libc.so.6 => /lib64/libc.so.6 (0x00007f9a7a46f000)
	libbz2.so.1 => /lib64/libbz2.so.1 (0x00007f9a7a25f000)
	libicuuc.so.50 => /lib64/libicuuc.so.50 (0x00007f9a79ee6000)
	libicui18n.so.50 => /lib64/libicui18n.so.50 (0x00007f9a79ae7000)
	libicudata.so.50 => /lib64/libicudata.so.50 (0x00007f9a78514000)
	/lib64/ld-linux-x86-64.so.2 (0x00007f9a7ca92000)
	libdl.so.2 => /lib64/libdl.so.2 (0x00007f9a78310000)

