#!/bin/bash

L1_R1= /media/vinbdi/data/tienanh/NSAIDS_pilot/NSAIDS_0011_L001_ds.81df2cc1573849138d82f8230c815a67/
NSAIDS_0011_L001_ds.81df2cc1573849138d82f8230c815a67  NSAIDS_0011_L003_ds.111556611a124592b9417420e504a5cf
NSAIDS_0011_L002_ds.be566a5456834aef871b632a3f6b327b  NSAIDS_0011_L004_ds.baf99a87df4d48008907ae34c2d0e2b3



data_folder='/media/vinbdi/data/NMT/Practice/FASTQ'
out_folder='/media/vinbdi/data/NMT/Practice/GATKout'
ref_gen='/home/vinbdi/Desktop/ref38'
ref_fasta='/home/vinbdi/Desktop/ref38/hg38.fasta'
dbsnp='/home/vinbdi/Desktop/ref38/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf'
b_mill='/home/vinbdi/ref38/Mills_and_1000G_gold_standard.indels.hg38.vcf'
folder='/media/vinbdi/data/NMT/Practice/FASTQ/SRR359199'
#target='/Users/jacky/work/VBDI-VAR/RB1_padding_100.bed'
PICARD='/media/vinbdi/data/NMT/Practice/picard/picard.jar'
# VARDICT='/home/vinbdi/VarDict-master'
#for folder in $data_folder/*; do
    sample=`(echo $folder | cut -d'/' -f8)`
    # echo "###############################"
    # echo $sample
    # check if folder exist
    des_folder=$out_folder/$sample
    # echo $des_folder
    # echo "###############################"

    #if [[ ! -d $des_folder ]]
    #then
        
        mkdir $out_folder/$sample
        mkdir $out_folder/$sample/fastqc
        mkdir $out_folder/$sample/aln
        mkdir $out_folder/$sample/qualification
        mkdir $out_folder/$sample/vcf
        mkdir $out_folder/$sample/annotation
         
        file1=`find $folder -type f \( -iname \*_1.filt.fastq.gz \)` #forward read -name is not case sensitive

        file2=`find $folder -type f \( -iname \*_2.filt.fastq.gz \)` #reverse read
        echo $file1
        echo $file2

        bwa mem -M \
            -t 20 \
            -R $(echo "@RG\tID:${sample}\tSM:${sample}\tLB:lib1\tPL:ILLUMINA") \
            $ref_fasta \
            $file1\
            $file2 | samtools view -Shb -o $out_folder/$sample/aln/${sample}_hg38_mapped.bam #b=output BAM, h=include header in output, S=Ignored for compatibility with previous samtools versions. Previously this option was required if input was in SAM format, but now the correct format is automatically detected by examining the first few characters of input.
        echo "@RG\tID:${sample}\tSM:${sample}\tLB:lib1\tPL:ILLUMINA"
        echo $out_folder/$sample/aln/${sample}_hg19_mapped.bam
        
        java -jar $PICARD FastqToSam \
            F1=$file1 \
            F2=$file2 \
            O=$out_folder/$sample/aln/${sample}_hg38_unmapped.bam \
            SAMPLE_NAME=`(echo $sample)` \
            READ_GROUP_NAME=`(echo $sample)` \
            LIBRARY_NAME=lib1 \
            PLATFORM=ILLUMINA
        #BAM stores more info, both forward and reverse reads

        java -jar $PICARD MergeBamAlignment \
            ALIGNED=$out_folder/$sample/aln/${sample}_hg38_mapped.bam \
            UNMAPPED=$out_folder/$sample/aln/${sample}_hg38_unmapped.bam \
            O=$out_folder/$sample/aln/${sample}_hg38_fin.bam \
            R=$ref_fasta
        #contain info from both aligned and unaligned bam files, can be used by Picard and GATK tools

        mark and sort duplicate file
        java -jar $PICARD MarkDuplicates I=$out_folder/$sample/aln/${sample}_hg38_fin.bam \
            O=$out_folder/$sample/qualification/${sample}_hg38_marked.bam \
            M=$out_folder/$sample/qualification/${sample}_hg38_metrics.txt \
            ASO=queryname \
            TMP_DIR=tmp
        
        java -jar $PICARD MarkDuplicates I=$out_folder/$sample/aln/${sample}_hg38_fin.bam \
                    O=$out_folder/$sample/qualification/${sample}_hg38_optical_marked.bam \
                    M=$out_folder/$sample/qualification/${sample}_hg38_optical_metrics.txt \
                    TAGGING_POLICY=OpticalOnly \
                    ASO=queryname \
                    TMP_DIR=tmp

        java -jar $PICARD MarkDuplicates I=$out_folder/$sample/aln/${sample}_hg38_fin.bam \
                    O=$out_folder/$sample/qualification/${sample}_hg38_marked_dupremoved.bam \
                    M=$out_folder/$sample/qualification/${sample}_hg38_metrics_dupremoved.txt \
                    REMOVE_DUPLICATES=true \
                    ASO=queryname \
                    TMP_DIR=tmp 
        #reads using REMOVE_DUPLICATES = #reads after filtering out Flag 0x0400 in the normal marked file because in this case there is no optical dup. In case there is optical dup, use also REMOVE_SEQUENCING_DUPLICATES to check

        sort sam: the index can only used after sorted based on coordinate
        java -jar $PICARD SortSam I=$out_folder/$sample/qualification/${sample}_hg38_marked.bam \
            O=$out_folder/$sample/qualification/${sample}_hg38_sorted.bam \
            SORT_ORDER=coordinate \
            TMP_DIR=tmp
        samtools index $out_folder/$sample/qualification/${sample}_hg38_sorted.bam

        # base calibration
        gatk BaseRecalibrator -I $out_folder/$sample/qualification/${sample}_hg38_sorted.bam \
            -R $ref_fasta --known-sites $dbsnp \
            -O $out_folder/$sample/qualification/${sample}_hg38_recal.table
        gatk ApplyBQSR -R $ref_fasta -I $out_folder/$sample/qualification/${sample}_hg38_sorted.bam \
            --bqsr-recal-file $out_folder/$sample/qualification/${sample}_hg38_recal.table \
            -O $out_folder/$sample/qualification/${sample}_hg38_arr.bam

        final_bam=$out_folder/$sample/qualification/${sample}_hg38_arr.bam

        # variant calling
        # HaplotypeCaller
        gatk HaplotypeCaller -R $ref_fasta \
            -I $final_bam \
            -O $out_folder/$sample/vcf/${sample}_HL_hg38.g.vcf \
            -ERC GVCF
        #-L $target \ after -R

        # gatk GenotypeGVCFs \
            -R $ref_fasta \
            -V $out_folder/$sample/vcf/${sample}_HL_hg38.g.vcf \
            -O $out_folder/$sample/vcf/${sample}_HL_hg38.vcf

        # Freebayes
        freebayes -f $ref_fasta \
            $final_bam>$out_folder/$sample/vcf/${sample}_FB_hg38.vcf
        # -t $target \
        
        #VarDict (haven't tested)
        AF_THR="0.01" # minimum allele frequency
        vardict -G $ref_fasta \
            -f $AF_THR -N $sample \
            -b $final_bam \
            -c 1 -S 2 -E 3 -g 4\
            $target|teststrandbias.R | var2vcf_valid.pl \
            -N $sample \
            -f $AF_THR>$out_folder/$sample/vcf/${sample}_VD_hg38.vcf



        # rm $out_folder/$sample/aln/${sample}_L1_hg19_fin.bam
        # rm $out_folder/$sample/aln/${sample}_L2_hg19_fin.bam
        # rm $out_folder/$sample/aln/${sample}_L1_hg19_mapped.bam
        # rm $out_folder/$sample/aln/${sample}_L2_hg19_mapped.bam
        # rm $out_folder/$sample/aln/${sample}_L1_hg19_unmapped.bam
        # rm $out_folder/$sample/aln/${sample}_L2_hg19_unmapped.bam
        # rm $out_folder/$sample/qualification/${sample}_hg19_sorted.bam
        # rm $out_folder/$sample/qualification/${sample}_hg19_marked.bam
        # rm $out_folder/$sample/aln/${sample}_hg19_merged_fin.bam

    #fi
#done
