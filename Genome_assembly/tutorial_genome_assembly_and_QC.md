### 1. Running FastQC on raw illumina data
* FastQC is a program that can quickly scan your raw data to help figure out if there are adapters or low quality reads present. Create a job file to run FastQC on one of the fastq files here: ```/data/genomics/workshops/smsc_2024/rawdata/```
	+ **module**: ```bio/fastqc```
	+ **command**: ```fastqc <FILE.fastq> -o .```
	+ after your job finishes, find the results and download some of the images, e.g. ```per_base_quality.png``` to your local machine using ffsend (load ```tools/ffsend``` module) and then the command ```ffsend upload <FILE>```.


### 2. Trimming adapters with TrimGalore! 
* PacBio data will be error-corrected solely by the assembler, but Illumina data trimming and thinning are common.
* Most assemblers these days don't want you to trim/thin for quality before assembling, but trimming is important for downstream applications. TrimGalore will auto-detect what adapters are present and remove very low quality reads (quality score <20) by default.  
* Create a job file to trim adapters and very low quality reads for the Illumina data here: ```/data/genomics/workshops/smsc_2024/rawdata/```
	+ **command**: ```trim_galore --paired --retain_unpaired <FILE_1.fastq> <FILE_2.fastq>```  
	+ **module**: ```bio/trim_galore```
	+ You can then run FastQC again to see if anything has changed.

### 3. Run Genomescope

* Genomescope can be used to estimate genome size from short read data: 
	+ [Genomescope](http://genomescope.org/genomescope2.0/) 

* To run Genomescope, first you need to generate a Jellyfish histogram.

* You'll need two job files for Jellyfish, one to count the kmers and the second to generate a histogram to give to Genomescope: 
* Here is a copy of the Guam Rail Illumina data: ```/data/genomics/workshops/smsc_2024/rawdata/```
	+ Hint: don't copy these data to your own space - they are very big.

* First job file: kmer count:
	+ Module: ```bio/jellyfish```
	+ Commands: ```zcat /data/genomics/workshops/smsc_2024/rawdata/SRR25828983_1.fastq.gz /data/genomics/workshops/smsc_2024/rawdata/SRR25828983_2.fastq.gz | jellyfish count -C -m 21 -t $NSLOTS -s 800000000 /dev/fd/0 -o reads.jf```
	+ ```zcat``` = uncompress the .fastq.gz files and sends the output to standard output (stdout)
 	+  ```/dev/fd/0``` = This allows jellyfish to read from the standard input (stdin), which is receiving the uncompressed data from zcat.
  	+ ```Pipping (|)``` = Passes the uncompressed data from zcat directly to the jellyfish count command.
  	+ ```jellyfish count``` = Command to Count k-mers
  	+ ```-m``` = kmer length  
	+ ```-s``` = sets the hash size. This depends on your system memory.
	+ ```-C``` = "canonical kmers" don't change this 
	+ Hint: this job needs to run on the high memory queue. 

* This will take a while, so we can move on and then come back to it when it finishes.

* Second job file: histogram:
	+ Module: ```bio/jellyfish```
	+ Commands: ```jellyfish histo -t $NSLOTS reads.jf > reads.histo```

* Download the histogram to your computer (e.g. using ffsend again), and put it in the Genomescope webservice: [Genomescope](http://genomescope.org/genomescope2.0/)

* let's run the analysis the same analysis with the HiFi PacBio data.

### 4. Hifiasm Assembly

Hifiasm is a fast haplotype-resolved de novo assembler for PacBio HiFi reads. It can produce partially phased assemblies of great quality. Hifiasm can use trio (Parental) short reads data or Hi-C data to produce haplotype-resolved assemblies.

#### Assembly using only HiFi data

* Job file: hifi_only.job

  + **PE:** multi-thread
  + **Queue:** Long, himem 
  + **Number of CPUs:** 32
  + **RAMMemory:** 10G (10G per CPU, 300G total)
  + **Module:** `module load bio/hifiasm`
  + **Command:**
```hifiasm -o Guam_Rail_only.asm -t 32 /data/genomics/workshops/smsc_2024/rawdata/SRR27030659_1_pacbio.fastq```

##### Comand explanation:
```
-o: name and path of the output file in asm format
-t: sets the number of CPUs in use

SRR27030659_1_pacbio.fastq are our PacBio Guam Rail data. Input sequences should be FASTA 
or FASTQ format, uncompressed or compressed with gzip (.gz). The quality scores of reads 
in FASTQ are ignored by hifiasm. Hifiasm outputs assemblies in `GFA <https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md>`_ format.
```
 
* This job should complete in a few hours

#### Assembly using both HiFi data and Hi-C

* Job file: hifi_Hi_C.job

  + **PE:** multi-thread
  + **Queue:** Long, himem 
  + **Number of CPUs:** 32
  + **RAMMemory:** 10G (10G per CPU, 300G total)
  + **Module:** `module load bio/hifiasm`
  + **Command:**
  ```hifiasm -o Guam_Rail_hic.asm -t32 --h1 7AFA556_1.fastq.gz --h2 7AFA556_2.fastq.gz Aambiguus_duplex.fastq.gz```

* This job should complete in a few hours
* In this mode, each contig represents a haplotig, which means that comes from one parental haplotype only.
* Hifiasm doen not performed scaffolding. For this run scaffolder such as SALSA or yahs.

##### Command explanation:
```
-o: name and path of the output file in asm format
-t: sets the number of CPUs in use
--h1: path to read 1 of Hi_C data
--h2: path to read 2 of Hi_C data
cloud_leopard_hifi.fq.gz Input reads. Input sequences should be FASTA 
or FASTQ format, uncompressed or compressed with gzip (.gz). The quality scores of reads 
in FASTQ are ignored by hifiasm. Hifiasm outputs assemblies in `GFA <https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md>`_ format.
```

One last step. convert your primary and other assemblies into a fasta files

```
qrsh
cd /path/to/your/assembly_results/
awk '/^S/{print ">"$2;print $3}' test.bp.p_ctg.gfa > test.p_ctg.fa
```

### 5. Run the fasta metadata parser to get statistics about both the primary assembly and haplotype-resolved assemblies.
* We use a python script to grab some statistics from assembly files. These are the stats it produces:  

Total number of base pairs:    
Total number of contigs:   
N90:  
N80:  
N70:  
N60:  
N50:  
L90:  
L80:  
L70:  
L60:  
L50:  
GC content:  
Median contig size:  
Mean contig size:  
Longest contig is:  
Shortest contig is: 

* We have put a finished assembly here: ```/data/genomics/workshops/smsc_2024/Guam_rail_assembly/bHypOws1_hifiasm.bp.p_ctg.fasta.gz```
	+ Module: ```bio/assembly_stats```
	+ Commands: ```assembly_stats <ASSEMBLY> > assembly_stats.out```

* How long is the longest contig and scaffold?
* What is the contig and scaffold N50?
* Are there differences between all the assemblies?
