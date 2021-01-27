#single sample
#from fastq

import os

out_dir = "/home/ted/ubuntu/ADR/out"
#os.mkdir(out_dir)

ref_fasta = 'path_ref_fasta'
dbsnp = 'path_dbsnp'
resource = 'path_to_resource'

G1_R1 = 'path11'
G1_R2 = 'path12'
G2_R1 = 'path21'
G2_R2 = 'path22'
G3_R1 = 'path31'
G3_R2 = 'path32'
G4_R1 = 'path41'
G4_R2 = 'path42'

sample = {}
sample['name']= "NSAIDS_0011"

sample['groups']=[]

sample['groups'].append({'groupname':'L001','read1':G1_R1,'read2':G1_R2})
sample['groups'].append({'groupname':'L002','read1':G2_R1,'read2':G2_R2})
sample['groups'].append({'groupname':'L003','read1':G3_R1,'read2':G3_R2})
sample['groups'].append({'groupname':'L004','read1':G4_R1,'read2':G4_R2})

out_aln = os.path.join(out_dir,'aln')
#os.mkdir(out_aln)


for group in sample['groups']:

    #map reads to reference
    group['mbam']= os.path.join(out_aln,sample['name']+'_'+group['groupname']+'_mapped.bam')
    cmd = f"""bwa mem -M \
                -t 20 \
                -R @RG\tID:${group['groupname']}\tSM:${sample['name']}\tLB:lib1\tPL:ILLUMINA \
                {ref_fasta} \
                {group['read1']} \
                {group['read2']} \
                | samtools view -Shb -o {group['mbam']}"""
    #os.system(cmd)

#merge, mark duplicate, sort

out_qual = os.path.join(out_dir,'qualification')
#os.mkdir(out_qual)

tmp = os.path.join(out_dir,'tmp')
#os.mkdir(tmp)

group['markedbam']= os.path.join(out_qual,sample['name']+'_marked.bam')
group['metrics']= os.path.join(out_qual,sample['name']+'_metrics.txt')

cmd = f"""gatk MarkDuplicates \
            -I {sample['groups'][0]['mbam']} \
            -I {sample['groups'][1]['mbam']} \
            -I {sample['groups'][2]['mbam']} \
            -I {sample['groups'][3]['mbam']} \
            -O {group['markedbam']} \
            -M {group['metrics']} \
            --TMP_DIR {tmp}"""
#os.system(cmd)

#base calibration
group['recaltable']= os.path.join(out_qual,sample['name']+'_recal.table')
cmd = f"""gatk BaseRecalibrator \
        -I {group['markedbam']} \
        -R {ref_fasta} \
        --known-sites {dbsnp} \
        -O {group['recaltable']}"""
#os.system(cmd)

group['arrbam']= os.path.join(out_qual,sample['name']+'_arr.bam')
cmd = f"""gatk ApplyBQSR \
            -R {ref_fasta} \
            -I {group['markedbam']} \
            --bqsr-recal-file {group['recaltable']} \
            -O {group['arrbam']}"""
#os.system(cmd)

#variant calling
out_vcf = os.path.join(out_dir,'vcf')
#os.mkdir(out_vcf)

group['vcf']= os.path.join(out_vcf,sample['name']+'.vcf')
group['bamout']= os.path.join(out_vcf,sample['name']+'out.bam')
cmd = f"""gatk HaplotypeCaller \
            -R {ref_fasta} \
            -I {group['arrbam']} \
            -O {group['vcf']} \
            -bamout {group['bamout']}"""
#os.system(cmd)

#CNN_filtering
out_fil = os.path.join(out_dir,'filter')
#os.mkdir(out_fil)

#1D_CNN_filter
group['1d_cnn_scored_vcf']= os.path.join(out_fil,sample['name']+'1d_cnn_scored.vcf')
cmd = f"""gatk CNNScoreVariants \
        -R {ref_fasta} \
        -V {group['vcf']} \
        -O {group['1d_cnn_scored_vcf']}"""
#os.system(cmd)

group['1d_cnn_filtered_vcf']= os.path.join(out_fil,sample['name']+'1d_cnn_filtered.vcf')
cmd = f"""gatk FilterVariantTranches \
        -V {group['1d_cnn_scored_vcf']} \
        -O {group['1d_cnn_filtered_vcf']} \
        --resource {resource} \
        --info-key CNN_1D \
        --snp-tranche 99.9 \
        --indel-tranche 95.0"""
#os.system(cmd)

#2D_CNN_filter
group['2d_cnn_scored_vcf']= os.path.join(out_fil,sample['name']+'2d_cnn_scored.vcf')
cmd = f"""gatk CNNScoreVariants \
        -R {ref_fasta} \
        -I {group['bamout']} \
        -V {group['vcf']} \
        -O {group['2d_cnn_scored_vcf']}
        --tensor-type read_tensor \
        --transfer-batch-size 8 \
        --inference-batch-size 8"""
#os.system(cmd)

group['2d_cnn_filtered_vcf']= os.path.join(out_fil,sample['name']+'2d_cnn_filtered.vcf')
cmd = f"""gatk FilterVariantTranches \
        -V {group['2d_cnn_scored_vcf']} \
        -O {group['2d_cnn_filtered_vcf']} \
        --resource {resource} \
        --info-key CNN_2D \
        --snp-tranche 99.9 \
        --indel-tranche 95.0"""
#os.system(cmd)