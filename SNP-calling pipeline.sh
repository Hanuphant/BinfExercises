#!/bin/bash

# Please install the following required dependencies before hand you use the script
# 1. BWA 
# 2. SAMTOOLS

# Before using the script I'd suggest you keep the reads 1, reads 2 file and reference genome file ref.fa in the directory where you will run the script

# Taking in your working directory for future use
directorypath=$(pwd)
##### 0.0 DECLARATION OF VARIABLES #####
# Initial declaration of the variables you wanted
## Think of passing it without generating error
#ref="~/class/ex6/hg38-chr17-refdir/"
#millsFile="~/class/ex6/millsFile/"
reads1=0
reads2=0
realign=0
gunzip=0
v=0
index=0
answer=0
outputstatus=0
mill=0
out=0

##### 0.1 USE OF GET-OPTS TO PASS ARGUMENTS #####
# Using getopts to pass arguments for this SNP pipeline
while getopts "a:b:r:eof:zvih" option; do
    case $option in
        a)  reads1="$OPTARG";;
        b)  reads2="$OPTARG";;
        r)  ref="$OPTARG";;
        e)  echo "Performing read re-alignment"; realign=1;index=1;;
        o)  out="$OPTARG";echo "Presenting an output file"; outputstatus=1;;
        f)  millsFile="$OPTARG";mill=1;;
        z)  echo "Compressing it to gunzip format"; gunzip=1;;
        v)  set -x; v=1;;
        i)  index=1;;
        h)  echo -e "\n---------Displaying Help---------\n"
            echo -e "This script is used to map FASTQ reads to the reference genome, realign and call the variants. For this package the following are the required dependencies: \n\t1. bwa for alignment \n\t2. samtools and bcftools for processing of the alignment map files. \nThis script will also use the GATK 3.7.0 for improving the alignment. Please note you don't have to install GATK. This script itself will install it.\nPass the following arguments to use the script.\n\t -a : Enter the complete path of reads1 \n\t -b : Enter the complete path of reads2 \n\t -r : Enter the complete path of reference genome\n\nPlease note that these arguments are necessary and my script won't work properly without them.\n\n\n\t -o : Generate an output of the script into a file \n\t -e : Perform realignment of the BAM file using GATK 3.7.0.\n\t -f : There is also a default Mills file taken from golden path. You can enter yours if you want to.\n\t -z : Generate and compress the output of my script in to compessed file of gunzip format.\n\t -v : Will run the script in verbose mode \n\t -i : Generate an indexing of the BAM file \n\t -h : You are already here. Help section for my script. \n\nMy script will also do file check but won't prompt you to enter the file again. It will just exit and ask you to try again. \n\nIf you still come across any issues please contact sgupta755@gatech.edu.\nScript is free to use and share.\n------------------------\n"; exit;;
    esac
done

### Initiating assignment

###### 0. FILE CHECKS ######

##### 0.1 CHECKS #####
if [ $reads1 -eq 0 ]
then
    echo "You haven't entered reads1. Please enter that."
    exit
fi
if [ $reads2 -eq 0 ]
then
    echo "You haven't entered reads2. Please enter that."
    exit
fi
if [ $ref -eq 0 ]
then
    echo "You haven't entered reference genome. Please enter that."
    exit
fi
if ! test -f "$reads1" 
then 
    echo "reads1 does not exist."
    exit
fi
if ! test -f "$reads2" 
then 
    echo "reads2 does not exist."
    exit
fi
if ! test -f "$ref" 
then 
    echo "Reference genome does not exist."
    exit
fi

###### 1. MAPPING OF READ FILE TO REF GENOME #######
#### 1.1 USING WGS WORKFLOW #####
# Telling bwa the path of reference genome

bwa index $ref

# Decompressing the gunzip input file 
gunzip -k -d -q $reads1
reads1=$(echo $reads1 | sed -E 's/.gz$//')

gunzip -k -d -q $reads2
reads2=$(echo $reads2 | sed -E 's/.gz$//')

#Making an output directory here
if ! test -d $directorypath/outputdirectory
then
    mkdir outputdirectory
fi

if test -f $directorypath/outputdirectory/output.vcf
then
    echo -e "There seems to be an output file created already. Do you wish to overwrite the output file (Y/n) ?"
    read fileprompt
    if [ $fileprompt == Y ]||[ $fileprompt == y ]
    then
        rm $directorypath/outputdirectory/output.vcf
    else 
        exit
    fi    
fi
#### 1.2 MAPPING VIA BWA ####
# Mapping via "bwa"
# This step is gonna take about 36 minutes for i7-11450H Quadcore - Single Thread
bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $ref $reads1 $reads2 -t 6 > $directorypath/outputdirectory/lane.sam 
echo -e "Mapping file saved in the $directorypath/outputdirectory/lane.sam.\n"
# read shreyash
# echo -e "Acknowledged"
#### 1.3 FIXING UP SAM FILE AND CONVERTING TO BAM FILE ####
samtools fixmate -O bam $directorypath/outputdirectory/lane.sam $directorypath/outputdirectory/lanefixed.bam -@ 6

# Makes a BAM file in 3-4 minutes
#### 1.4 CREATION OF FASTA INDEX AND SEQUENCE DICTIONARY FILE ####
# Creates FASTA index file 
samtools faidx $ref
# Creates sequence dictionary file
len=$(expr length $ref)
refpath=${ref:0:$len-4}
samtools dict $ref -o $refpath.dict

#$nameoftheReferenceFile=$(echo $ref|sed "s:$directorypath/::")
# samtools dict $ref -o $directorypath/$nameoftheReferenceFile.dict

if [ $index -eq 1 ]
then
    #### 1.5 SORTING THE BAM FILE AS GATK REQUIRES INDEXED BAM FILE AS INPUT ####
    samtools sort -O bam -o $directorypath/outputdirectory/lane_fixed_sorted.bam -T /tmp/lane_temp $directorypath/outputdirectory/lanefixed.bam -@ 6
    # Indexing the sorted BAM file
    samtools index -b $directorypath/outputdirectory/lane_fixed_sorted.bam -@ 6
fi

###### 2. REALIGNMENT USING GENOME ANALYSIS TK ######

if [ $realign -eq 1 ]
then
    ##### 2.1 DOWNLOAD AND DECOMPRESS #####
    echo "Need JAVA installed for further process."

    # 2.1.1 Downloading GATK 3.7.0 from Google Cloud
    wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2
    
    # Decompressing it
    tar -xjpf GenomeAnalysisTK-3.7-0-gcfedb67.tar.bz2

    ##### 2.2 DOWNLOADING MILLS FILE FROM THE DATABASE SERVER #####
    if [ $mill -eq 0 ]
    then
        # 2.2.1 Download command
        wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        # 2.2.2 Extracting it via gunzip
        gzip -d -k -q Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
        millsFile="$directorypath/Mills_and_1000G_gold_standard.indels.hg38.vcf"
    fi

    ##### 2.3 REALIGNMENT FROM GENOME ANALYSIS TK 3.7.0 #####

    # Creating target intervals for the output BAM file
    java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I $directorypath/outputdirectory/lane_fixed_sorted.bam -o $directorypath/outputdirectory/lane.intervals --known $millsFile 2> sgupta755.log
    # Took about to 132.77 secs

    # Using target interval creating realigned BAM file
    java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I $directorypath/outputdirectory/lane_fixed_sorted.bam -targetIntervals $directorypath/outputdirectory/lane.intervals -known $millsFile -o $directorypath/outputdirectory/lane_fixed_sorted_realigned.bam 2>> sgupta755.log
    # Took about 192.43 secs

    # After indexing now is the time for index our BAM file again
    samtools index -b $directorypath/outputdirectory/lane_fixed_sorted_realigned.bam

else
    if [ index -eq 1 ]
    then
        cp $directorypath/outputdirectory/lane_fixed_sorted.bam $directorypath/outputdirectory/lane_fixed_sorted_realigned.bam
    fi 
    cp $directorypath/outputdirectory/lanefixed.bam $directorypath/outputdirectory/lane_fixed_sorted_realigned.bam
fi

###### 3. VARIANT CALLING ######
# Remember from the alignment we generated one BAM file
##### 3.1 CREATING A VCF FILE #####
bcftools mpileup -A -Ou -f $ref $directorypath/outputdirectory/lane_fixed_sorted_realigned.bam|bcftools call -vmO z -o $directorypath/outputdirectory/lane_fixed_sortedvcf.vcf.gz
# This takes a bit long with no progess bar
# Okay very long

##### 3.2 INDEXING THE VCF FILE USING BCFTOOLS #####
bcftools index $directorypath/outputdirectory/lane_fixed_sortedvcf.vcf.gz

##### 3.3 DIDN'T DO FILTERING OF VARIANTS #####


###### 4. CREATING BED FILE FROM VCF FILE ######
##### 4.1 REMOVING HEADERS #####
# Removing headers via sed by substituting anything which starts with # to nothing
gzip -d -q $directorypath/outputdirectory/lane_fixed_sortedvcf.vcf.gz
cd $directorypath/outputdirectory/
sed -ie "/^##/d" lane_fixed_sortedvcf.vcf

# Seperating columns from 1-5
cut -f1-5 lane_fixed_sortedvcf.vcf > edited_cut_lane_fixed_sorted_vcf.vcf
# Creating values for LENGTH 
awk 'BEGIN{FS=OFS="\t"}{x=$4;y=$5;$6=length(y)-length(x);print $0}' edited_cut_lane_fixed_sorted_vcf.vcf > edit1.vcf
# Creating values for STOP
awk 'BEGIN{FS=OFS="\t"}{x=$2;y=$6;$7=x+y;print $0}' edit1.vcf > edit2.vcf
# Selecting CHROM, POS, LENGTH, STOP
cut -f1,2,6,7 edit2.vcf > edit3.vcf
# Exchanging the rows
awk 'BEGIN{FS=OFS="\t"}{x=$3;y=$4;$3=y;$4=x;print $0}' edit3.vcf > edit4.vcf
# Renaming of columns
sed -ie "s:POS\t0\t0:START\tSTOP\tLENGTH:" edit4.vcf
# Removing chr from chr17
sed -ie "s:^chr::" edit4.vcf
awk 'BEGIN{FS=OFS="\t"}{if ($4 ~ "LENGTH" || $4==0){print $0 > ("snps.txt")}}{if ($4~"LENGTH" || $4!=0){print $0 > ("indels.txt")}}' edit4.vcf
# The variant is an SNP if there is no change bases
# The variant is an indel if there is a change in bases or rather there is an insertion or deletion.
cat edit4.vcf > output.vcf
if [ $outputstatus -eq 1 ]
then
    cp $directorypath/outputdirectory/lane_fixed_sortedvcf.vcf.gz $out
    gzip -d $out
    if [ $gunzip -eq 1 ]
    then 
        gzip $out
    fi
    rm output.vcf
fi
if [ $outputstatus -ne 1 ]
then
    rm output.vcf
fi
##### CLEANER #####
rm edit*

exit