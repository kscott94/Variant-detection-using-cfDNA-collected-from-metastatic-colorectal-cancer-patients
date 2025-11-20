#!/bin/bash

### Software ###################################################################
#
# SRA-Toolkit v3.2.1
# FastQC v0.11.9 
# MultiQC v1.32
# Samtools v1.21
# BCFtools v1.20
# BWA v0.7.17-r1188
# GATK v4.6.2.0
# VarScan v2.3.9
# Docker v26.1.3
# Ensembl-VEP v115.2
# Jupyter-notebook v
# R v4.5.1
#



### Set up #####################################################################
proj_dir="/compute_space/GT_project/submission"
cd ${proj_dir}

proj_dir="/contents"
#Install software
wget https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip
wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar/download -O VarScan.v2.3.9.jar

alias gatk="/source_files/gatk-4.6.2.0/gatk"
alias varscan2="java -jar /source_files/VarScan.v2.3.9.jar"

#Pull VEP docker image (variant annotation tool)
docker pull ensemblorg/ensembl-vep


#Gather reference material

#Get reference genome from broad institute and generate indices
mkdir ${proj_dir}/references
cd ${proj_dir}/references

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

#fasta index
samtools faidx ${proj_dir}/references/Homo_sapiens_assembly38.fasta

#BWA index (This takes like 3 hours)
bwa index ${proj_dir}/references/Homo_sapiens_assembly38.fasta

#GATK dictionary
gatk CreateSequenceDictionary -R ${proj_dir}/references/Homo_sapiens_assembly38.fasta

#Get panel of normals from 1000 genomes project and other natural variants
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz -O dbsnp_155.hg38.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz.tbi -O dbsnp_155.hg38.vcf.gz.tbi
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi


## Get fastq files for patient CTDC42 (CTC276)
mkdir ${proj_dir}/reads
cd ${proj_dir}/reads

fasterq-dump SRR13974004 --split-files #baseline
fasterq-dump SRR13974005 --split-files #first response
fasterq-dump SRR13974006 --split-files #PD

#Notes:
#patient IDs reported in supplementary files are not consisent with SRA library 
#IDs or sample names
#cross-matched number of sequenced bases and number of mapped reads between fastq 
#files and what was reported for each patient ID in supp table 2 to figure out 
#which SRA files belong to patient CTDC42
#raw reads were not deposited into SRA, only post-QC, mapped reads available

# Generate fastq QC reports of fastq files
fastqc *
multiqc .

#Notes:
#reads are already trimmed and QC'ed but still contain PCR duplicates?
#reads were filtered such that some reads are unpaired



### Map reads ##################################################################


## Map reads to human genome assembly 38 and remove PCR duplicates (authors used hg19)
#I will also run the pipeline on reads without deduplication since authors did 
#not remove PCR duplicates (to my knowledge).

#map reads using BWA and deduplicate
for R1 in *_1.fastq.gz
do
    R2=${R1/_1/_2} #get file name of mate 2
    SAMPLE=$(basename ${R1} _1.fastq.gz)  #get sample name by removing mate id and file extension
    REF=${proj_dir}/references/Homo_sapiens_assembly38.fasta
    
    echo "Processing ${SAMPLE}"
    
    bwa mem ${REF} ${R1} ${R2} \
    | samtools sort -n - \
    | samtools fixmate -m - - \
    | samtools sort - \
    | samtools markdup -r - ${SAMPLE}.rmdup.bam #detect and remove PCR duplicates
    
    samtools index ${SAMPLE}.rmdup.bam
   
done

#add sample name to bam read group 
for bam in *.rmdup.bam
do
    SAMPLE=$(basename ${bam} .rmdup.bam)  #get sample name by removing mate id and file extension
    gatk AddOrReplaceReadGroups \
        -I $bam \
        -O ${SAMPLE}.rmdup.rg.bam \
        -RGID ${SAMPLE} \
        -RGSM ${SAMPLE} \
        -RGLB lib1 \
        -RGPL ILLUMINA \
        -RGPU unit1

    samtools sort -o ${SAMPLE}.rmdup.rg.csorted.bam ${SAMPLE}.rmdup.rg.bam
    samtools index ${SAMPLE}.rmdup.rg.csorted.bam
    
    echo "################  DONE: ${SAMPLE} ##################"
done

#base recalibration
for bam in *.rmdup.rg.csorted.bam
do
    SAMPLE=$(basename $bam .rmdup.rg.csorted.bam)
    OUT_File=${SAMPLE}.rmdup.rg.csorted.bsrq.bam
    
    REF=${proj_dir}/references/Homo_sapiens_assembly38.fasta
    gnomAD=${proj_dir}/references/af-only-gnomad.hg38.vcf.gz #known/natural variants
    MILLS=${proj_dir}/references/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz #known/natural variants
    
    echo "Running BQSR for ${SAMPLE}"

    gatk BaseRecalibrator \
        -I ${bam} \
        -R ${REF} \
        --known-sites ${gnomAD} \
        --known-sites ${MILLS} \
        -O ${SAMPLE}.recal_data.table


    echo "Applying BQSR for ${SAMPLE}"

    gatk ApplyBQSR \
        -R ${REF} \
        -I ${bam} \
        --bqsr-recal-file ${SAMPLE}.recal_data.table \
        -O $OUT_File

    samtools index $OUT_File

    echo "############ DONE: ${SAMPLE} ############"
done


#generate BAM QC reports and verify PCR duplicates were removed
fastqc *
multiqc .



### Call Variants ##############################################################

# Due to limited computing power and time, only look at chromosome 17
# focus on somatic short mutations (SNVs/indels) only 
# no matched normal sample (cant find PBMC sample for this patient??), technical artifacts may be difficult to spot unless I analyze additional samples from other patients
# Since this patient does not have an assocated PBMC sample, I will have to filter out suspected germline variants manually (VAF<40%?)
# will use panel of normals from 1000 genomes project
# non-FFPE stample, so no need to model read orientation biase

## GATK pipeline -- Variant detection

for bam in *.rmdup.rg.csorted.bsrq.bam
do
    REF=${proj_dir}/references/Homo_sapiens_assembly38.fasta
    PON=${proj_dir}/references/1000g_pon.hg38.vcf.gz
    SAMPLE=$(basename ${bam} .rmdup.rg.csorted.bsrq.bam)
    OUT=${SAMPLE}.dedup.mutec2.chr17

    echo "Running Mutect2 for ${SAMPLE}"

    gatk Mutect2 \
        -R ${REF} \
        -I ${bam} \
        -tumor ${SAMPLE} \
        --panel-of-normals ${PON} \
        -L chr17 \
        -O ${OUT}.unfiltered.vcf.gz


    echo "Filtering low quality calls"
    
    gatk FilterMutectCalls \
        -V ${OUT}.unfiltered.vcf.gz \
        -R ${REF} \
        -O ${OUT}.filtered.vcf.gz


    echo "############# DONE: ${SAMPLE} ###############"
done



## VarScan2 pipeline -- Variant detection
#authors used Varscan2 for variant detection. So I will also.
#Min coverage:	8
#Min var freq:	0.01
#Min avg qual:	15
#P-value thresh:	0.05


for bam in *.rmdup.rg.csorted.bsrq.bam
do
    REF=${proj_dir}/references/Homo_sapiens_assembly38.fasta
    SAMPLE=$(basename ${bam} .rmdup.rg.csorted.bsrq.bam)
    OUT=${SAMPLE}.dedup.varscan.chr1
    
    samtools mpileup -r chr17 -f ${REF} ${bam} > ${SAMPLE}.mpileup
    
    varscan2 mpileup2snp ${SAMPLE}.mpileup \
        --output-vcf 1 --min-var-freq 0.01 \
        --p-value 0.05 > ${OUT}.snp.vcf
            
    varscan2 mpileup2indel ${SAMPLE}.mpileup \
        --output-vcf 1 --min-var-freq 0.01 \
        --p-value 0.05 > ${OUT}.indel.vcf

done

mv *vcf.gz ../../../vcf
mv *vcf ../../../vcf


### Annotate Variants ##########################################################

#annotate variants using

#move vcf files to new directory
proj_dir="/compute_space/GT_project/submission"
REF=${proj_dir}/references/Homo_sapiens_assembly38.fasta
mkdir ${proj_dir}/vcf
mkdir ${proj_dir}/vcf/annot
cd ${vcf_dir}/annot


#get functional annotation resource
gatk FuncotatorDataSourceDownloader \
  --somatic \
  --validate-integrity \
  --extract-after-download \
  --hg38
  -O ${REF}

funcotator_data_sources=${proj_dir}/references/


#annotate VCFs
for vcf in ../*.vcf.gz
do  
    SAMPLE=$(basename ${vcf} .vcf.gz)
    REF=${proj_dir}/references/Homo_sapiens_assembly38.fasta
    OUT=${SAMPLE}.annot
    
    gatk IndexFeatureFile -I ${vcf}
    
    echo "Annotating ${SAMPLE}"
    
    gatk Funcotator \
        -V ${vcf} \
        -R ${REF} \
        --ref-version hg38 \
        --data-sources-path ${funcotator_data_sources}/ \
        -O ${OUT}.vcf.gz \
        --output-file-format VCF
done   


 
    
    

