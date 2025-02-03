# ROHs

Runs of Homozygosity are stretches of DNA where an individual inherits identical haplotypes from both parents. They often are a result of inbreeding, which can occur when closely related individuals mate. 

We will use BCFtools roh to estimate Runs of Homozygosity


## Setting us up for analyses

First, lets make a new directory for our analyses. 

```
mkdir /scratch/genomics/YOURNAME/smsc_2024/geneticdiversity/
mkdir /home/YOURNAME/smsc_2024/day6
```

Next, we will copy all of the already written scripts into this folder so we can slightly alter it: 
```
cp /data/genomics/workshops/smsc_2024/scripts_Day6/*.job /home/YOURNAME/smsc_2024/day6
```


## Check depth of the bam files 

First, we need to check the average depth of our bam files. This is important, as to do an ROH analysis wth BCFtools, it requies genomes to be at least ~15x to confidently assess runs of homozygosity. 


Fulica_atra


```
# /bin/sh                                                                                                                                                               
# ----------------Parameters---------------------- #                                                                                                                    
#$ -S /bin/sh                                                                                                                                                           
#$ -pe mthread 5                                                                                                                                                        
#$ -q lThM.q                                                                                                                                                            
#$ -l mres=10G,h_data=5G,h_vmem=5G,himem                                                                                                                                
#$ -cwd                                                                                                                                                                 
#$ -j y                                                                                                                                                                 
#$ -N checkdepth                                                                                                                                                        
#$ -o checkdepth.log                                                                                                                                                    
#                                                                                                                                                                       
# ----------------Modules------------------------- #                                                                                                                    
module unload gcc/8.5.0                                                                                                                                                 
module load bio/samtools                                                                                                                                                
#                                                                                                                                                                       
# ----------------Your Commands------------------- #                                                                                                                    
#                                                                                                                                                                       
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                           
echo + NSLOTS = $NSLOTS                                                                                                                                                 
#                                                                                                                                                                       
#                                                                                                                                                                       
BAM=/data/genomics/workshops/smsc_2024/BAMS/HOW_N23-0568.realigned.bam                                                 

OUT=HOW_N23-0568_depth                                                                                           
#                                                                                                 
samtools depth ${BAM}  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${OUT}

echo = `date` job $JOB_NAME done                                                                  
#                                                                                                  
#                                                                                           
```
This will take a couple minutes to run. 

We can look at our per sample depth here: 
```
HOW_N23-0063 25.16
HOW_N23-0068 22.2183
Atlantisia_rogersi 19.52
Fulica_atra 37.33
Gallinula_chloropus 22.40
Grus_americana 27.45
Heliornis_fulica 4.43
Hypotaenidia_okinawae 70.48 
Laterallus_jamaicensis_coturniculus 11.28
Porphyio_hochstetteri 12.81
Rallus_limicola 15.01
Zapornia_atra 0.14


```

## Subset and filter by individual 

Next we will subset and filter each individual by its specific depth. 

We alse only include SNPs that had a depth of more than one third and less than double of the average depth of coverage for each sample

So we will only keep samples that are over 15x: 
```
HOW_N23-0063 25.16 
HOW_N23-0068 22.2183
Atlantisia_rogersi 19.52
Rallus_limicola 15.01
```

Now we can run the code: 
```
# /bin/sh                                                                                                                                                               
# ----------------Parameters---------------------- #                                                                                                                    
#$ -S /bin/sh                                                                                                                                                           
#$ -pe mthread 5                                                                                                                                                        
#$ -q lThM.q                                                                                                                                                            
#$ -l mres=10G,h_data=5G,h_vmem=5G,himem                                                                                                                                
#$ -cwd                                                                                                                                                                 
#$ -j y                                                                                                                                                                 
#$ -N filterbyind                                                                                                                                                       
#$ -o filter_indv.log                                                                                                                                                  
#                                                                                                                                                                       
# ----------------Modules------------------------- #                                                                                                                    
module unload gcc/8.5.0                                                                                                                                                 
module load bio/vcftools                                                                                                                                                
#                                                                                                                                                                       
# ----------------Your Commands------------------- #                                                                                                                    
#                                                                                                                                                                       
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                           
echo + NSLOTS = $NSLOTS                                                                                                                                                 
#                                                                                                                                                                       
#                                                                                                                                                                       
VCF=/data/genomics/workshops/smsc_2024/VCF/gatk_allsamples_ptg000001l_filtered.recode.vcf                                                                               
OUTDIR=/scratch/genomics/YOURNAME/smsc_2024/geneticdiversity/                                                                                                          
#                                                                                                                                                                       
NAME=GuamRail_HOW_N23_0063                                                                                                                                                    
MAX=75.48                                                                                                                                                               
MIND=12.58                                                                                                                                                              
#                                                                                                                                                                       
#                                                                                                                                                                       
vcftools --vcf ${VCF} --minDP ${MIND} --maxDP ${MAX} --indv ${NAME} \                                                                                                   
--min-alleles 2 --max-alleles 2 --max-missing 1 \                                                                                                                       
--out ${OUTDIR}/gatk_allsamples_ptg000001l_filtered_min${MIND}_max${MAX}_${NAME} \                                                                                      
--recode --recode-INFO-all                                                                                                                                              
#                                                                                                                                                                       
echo = `date` job $JOB_NAME done  
```

Let's look at the log file: 

```
less filter_indv.log
```
And then the output file: 
```
less -S /data/genomics/workshops/smsc_2024/VCF/gatk_allsamples_ptg000001l_filtered_min12.58_max75.48_GuamRail_HOW_N23_0063.recode.vcf
```


## Run BCFtools!

BCFtools roh is a program that estimates ROH.  We will use it now to estimate ROHs in our samples. 

```
# /bin/sh                                                                                                                                                               
# ----------------Parameters---------------------- #                                                                                                                    
#$ -S /bin/sh                                                                                                                                                           
#$ -pe mthread 5                                                                                                                                                        
#$ -q lThM.q                                                                                                                                                            
#$ -l mres=10G,h_data=5G,h_vmem=5G,himem                                                                                                                                
#$ -cwd                                                                                                                                                                 
#$ -j y                                                                                                                                                                 
#$ -N bcftools                                                                                                                                                          
#$ -o bcftools.log                                                                                                                                                      
#                                                                                                                                                                       
# ----------------Modules------------------------- #                                                                                                                    
module load bio/bcftools/1.19                                                                                                                                           
#                                                                                                                                                                       
# ----------------Your Commands------------------- #                                                                                                                    
#                                                                                                                                                                       
echo + `date` job $JOB_NAME started in $QUEUE with jobID=$JOB_ID on $HOSTNAME                                                                                           
echo + NSLOTS = $NSLOTS                                                                                                                                                 
#                                                                                                                                                                       
#                                                                                                                                                                       
VCF=/data/genomics/workshops/smsc_2024/VCF/gatk_allsamples_ptg000001l_filtered_min12.58_max75.48_GuamRail_HOW_N23_0063.recode.vcf                                                                                                                                                                       
#                                                                                                                                                                       
bcftools roh -G30 --AF-dflt 0.4 ${VCF} -o GuamRail_HOW_N23_0063.txt            
#                                                                                                                                                                       
echo = `date` job $JOB_NAME done  
```

Here are some meanings to the filters: 

- bcftools roh: This command runs the roh plugin from bcftools to detect runs of homozygosity.
- AF-dflt 0.4: This option sets the default allele frequency to 0.4 when the frequency is missing in the input file.
- G30: This option sets the phred-scaled genotype quality threshold to 30. Genotypes below this quality threshold will be treated as missing.


Let's look at the output file, it should run fast. 

```
less GuamRail_HOW_N23_0063.txt 
```

We can see regions where there are no ROH, and then when there regions where there are with RG.  

## Pull RG Flag 

Next we will pull out the RG flag from our samples. 

```
grep "RG" GuamRail_HOW_N23_0063.txt > GuamRail_HOW_N23_0063_RG.txt
```

Let's just look at the number of ROHs in each sample we have here. 

```
less GuamRail_HOW_N23_0063_RG.txt
```


## Plot ROH in R
```

library (ggplot2)
library(tidyverse)
library(dplyr)

############################################
### Plotting bar graph of ROH categories ###
############################################

#read in data
dat_0.1 <- read.csv ("RunsofHomozygosity.csv", header=TRUE)
nrow(dat_0.1)
#7338
##Remove ROH under 100kb (100kb=100,000bp)
dat <- subset(dat, Length.bp. > 100000)
nrow(dat)
#472

##################################################
##################################################
#Let's look at the number of ROH in each sample ## 
##################################################
##################################################
dat_Atlantisia_rogersi <- subset(dat, Sample == "Atlantisia_rogersi")
nrow(dat_Atlantisia_rogersi) 
#183
dat_GuamRail_HOW_N23_0063 <- subset(dat, Sample == "GuamRail_HOW_N23_0063")
nrow(dat_GuamRail_HOW_N23_0063) 
#61
dat_GuamRail_HOW_N23_0068 <- subset(dat, Sample == "GuamRail_HOW_N23_0068")
nrow(dat_GuamRail_HOW_N23_0068) 
#127
dat_Rallus_limicola <- subset(dat, Sample == "Rallus_limicola")
nrow(dat_Rallus_limicola) 
#101

#Let's make a bar graph for length vs. number of ROH

ROH_lengthvsnumber <- read.csv("ROH_lengthvsNumber.csv", header=TRUE)

 ggplot(ROH_lengthvsnumber,aes(x=Total_length_ROH,y=NumberofROH,label=Sample, color=Sample)) + geom_point(size = 3) + theme_classic()

#What different patterns do we see? 

######################################################
######################################################
### Let's make a graph of different lengths of ROH ###
######################################################
######################################################

##################################
# Get ROH values between 0.1-1Mb
##################################
dat_0.1_1Mb <- subset(dat, Length.bp. < 1000000)
sum <- aggregate(dat_0.1_1Mb$Length.bp., by=list(Category=dat_0.1_1Mb$Sample), FUN=sum)
sum

#               Category        x
#1    Atlantisia_rogersi 42456648
#2 GuamRail_HOW_N23_0063 13222867
#3 GuamRail_HOW_N23_0068 28095326
#4       Rallus_limicola 18974775

write.csv (dat_0.1_1Mb , "samples_between_0.1_1mb.csv")
##################################
# Get ROH values between 1-5Mb
##################################
dat_1_10Mbtest <- subset(dat, Length.bp. > 1000000)
dat_1_10Mb <- subset(dat_1_10Mbtest, Length.bp. < 5000000)
sum_2_5Mb <- aggregate(dat_1_10Mb$Length.bp., by=list(Category=dat_1_10Mb$Sample), FUN=sum)
sum_2_5Mb

#               Category       x
#1 GuamRail_HOW_N23_0063 5020811
#2       Rallus_limicola 1055016


##################################
# Get ROH values between 5-10Mb
##################################
dat_10_100Mb <- subset(dat, Length.bp. > 5000000)
dat_10_100Mbfinal <- subset(dat_10_100Mb, Length.bp. < 10000000)
sum_5_10Mb <- aggregate(dat_10_100Mbfinal$Length.bp., by=list(Category=dat_10_100Mbfinal$Sample), FUN=sum)
sum_5_10Mb

#               Category        x
#1 GuamRail_HOW_N23_0063 20835957

##################################
# Get ROH values between 10-100Mb
##################################
dat_10_100Mb <- subset(dat, Length.bp. > 10000000)
dat_10_100Mbfinal <- subset(dat_10_100Mb, Length.bp. < 100000000)
sum_10_100Mb <- aggregate(dat_10_100Mbfinal$Length.bp., by=list(Category=dat_10_100Mbfinal$Sample), FUN=sum)

#no rows 

######################
#Plotting by category #
######################

datrohlength <-read.csv("ROH_byLength_Category.csv", header=TRUE)

p <- ggplot(datrohlength, aes(fill=ROHcat, y=Length, x=reorder(Sample,-Totalsize ))) + 
    geom_bar(position="stack", stat="identity", color="black") + theme_classic() + coord_flip() + scale_fill_manual(values = c( "#e6ab02", "#d95f02" ,"#7570b3")) + scale_x_discrete(labels=c("5.0e+08" = "0.5", "1.0e+09" = "1","1.5e+09" = "1.5"))+  theme(axis.text.x = element_text(color="black", size = 10), axis.text.y = element_text(color="black", size = 10))
p




######################################################
######################################################
############# Plot cummulative ROH ###################
######################################################
######################################################

dat <- read.csv("RunsofHomozygosity_forcummulative_final.csv", header=TRUE)

#Filter out ROH
dat <- subset(dat, Length.bp. > 100000)
dat <- subset(dat, Quality > 80)

#Plotting
p <- ggplot(dat, aes(x=Length.bp., y=Cummulative, color=Sample)) + theme_classic() + scale_x_log10()   +
  geom_line(size=1.2, alpha=0.75) + scale_color_manual(values=c("orangered2", "mediumseagreen", "royalblue", "orange", "mediumseagreen" ))
p
q <- p + geom_line(aes(size = Bold))  +
  scale_size_manual(values = c(0.1, 1.5))  
q 


```
