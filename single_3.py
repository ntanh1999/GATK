#single sample
#Data processing + Calling variants

import os

#make directory 
out_dir = "/media/vinbdi/data/tienanh/gatk"
#os.mkdir(out_dir)
out_aln = os.path.join(out_dir,'aln')
#os.mkdir(out_aln)
out_qual = os.path.join(out_dir,'qualification')
#os.mkdir(out_qual)
tmp = os.path.join(out_dir,'tmp')
#os.mkdir(tmp)
out_vcf = os.path.join(out_dir,'vcf')
#os.mkdir(out_vcf)
out_fil = os.path.join(out_dir,'filter')
#os.mkdir(out_fil)

#create variables
ref_fasta = '/home/vinbdi/Desktop/ref38/hg38.fasta'
dbsnp = '/home/vinbdi/Desktop/ref38/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf'
resource = '???'

G1_R1 = '/media/vinbdi/data/tienanh/NSAIDS_pilot/NSAIDS_0011_L001_ds.81df2cc1573849138d82f8230c815a67/NSAIDS-0011_S1_L001_R1_001.fastq.gz'
G1_R2 = '/media/vinbdi/data/tienanh/NSAIDS_pilot/NSAIDS_0011_L001_ds.81df2cc1573849138d82f8230c815a67/NSAIDS-0011_S1_L001_R2_001.fastq.gz'
G2_R1 = '/media/vinbdi/data/tienanh/NSAIDS_pilot/NSAIDS_0011_L002_ds.be566a5456834aef871b632a3f6b327b/NSAIDS-0011_S1_L002_R1_001.fastq.gz'
G2_R2 = '/media/vinbdi/data/tienanh/NSAIDS_pilot/NSAIDS_0011_L002_ds.be566a5456834aef871b632a3f6b327b/NSAIDS-0011_S1_L002_R2_001.fastq.gz'
G3_R1 = '/media/vinbdi/data/tienanh/NSAIDS_pilot/NSAIDS_0011_L003_ds.111556611a124592b9417420e504a5cf/NSAIDS-0011_S1_L003_R1_001.fastq.gz'
G3_R2 = '/media/vinbdi/data/tienanh/NSAIDS_pilot/NSAIDS_0011_L003_ds.111556611a124592b9417420e504a5cf/NSAIDS-0011_S1_L003_R2_001.fastq.gz'
G4_R1 = '/media/vinbdi/data/tienanh/NSAIDS_pilot/NSAIDS_0011_L004_ds.baf99a87df4d48008907ae34c2d0e2b3/NSAIDS-0011_S1_L004_R1_001.fastq.gz'
G4_R2 = '/media/vinbdi/data/tienanh/NSAIDS_pilot/NSAIDS_0011_L004_ds.baf99a87df4d48008907ae34c2d0e2b3/NSAIDS-0011_S1_L004_R2_001.fastq.gz'

sample = {}
sample['name']= "NSAIDS_0011"
sample['groups']=[]
sample['groups'].append({'groupname':'L001','read1':G1_R1,'read2':G1_R2})
sample['groups'].append({'groupname':'L002','read1':G2_R1,'read2':G2_R2})
sample['groups'].append({'groupname':'L003','read1':G3_R1,'read2':G3_R2})
sample['groups'].append({'groupname':'L004','read1':G4_R1,'read2':G4_R2})


for group in sample['groups']:
    #map reads to reference
    group['mappedbam']= os.path.join(out_aln,sample['name']+'_'+group['groupname']+'_mapped.bam')
    cmd = f"""bwa mem -M \
                -t 20 \
                {ref_fasta} \
                {group['read1']} \
                {group['read2']} \
                | samtools view -Shb -o {group['mappedbam']}"""
    print(f'RUNNING {cmd}')
    #os.system(cmd)
    
    #add readgroup metadata
    group['addedbam']= os.path.join(out_aln,sample['name']+'_'+group['groupname']+'_added.bam')
    cmd = f"""gatk AddOrReplaceReadGroups \
            -I {group['mappedbam']} \
            -O {group['addedbam']} \
            --RGID {group['groupname']} \
            --RGSM {sample['name']} \
            --RGLB lib1 \
            --RGPU unit1 \
            --RGPL ILLUMINA"""
    print(f'RUNNING {cmd}')
    #os.system(cmd)
            
    #sortbam
    group['sortedbam']= os.path.join(out_aln,sample['name']+'_'+group['groupname']+'_sorted.bam')
    cmd = f"""gatk SortSam \
            -I {group['addedbam']} \
            -O {group['sortedbam']} \
            -SORT_ORDER coordinate \
            --TMP_DIR {tmp}"""
    print(f'RUNNING {cmd}')
    #os.system(cmd)

#mark duplicate
sample['markedbam']= os.path.join(out_qual,sample['name']+'_marked.bam')
sample['metrics']= os.path.join(out_qual,sample['name']+'_metrics.txt')

cmd = f"""gatk MarkDuplicates \
            -I {sample['groups'][0]['sortedbam']} \
            -I {sample['groups'][1]['sortedbam']} \
            -I {sample['groups'][2]['sortedbam']} \
            -I {sample['groups'][3]['sortedbam']} \
            -O {sample['markedbam']} \
            -M {sample['metrics']} \
            --TMP_DIR {tmp}"""
print(f'RUNNING {cmd}')
#os.system(cmd)

#base calibration
sample['recaltable']= os.path.join(out_qual,sample['name']+'_recal.table')
cmd = f"""gatk BaseRecalibrator \
        -I {sample['markedbam']} \
        -R {ref_fasta} \
        --known-sites {dbsnp} \
        -O {sample['recaltable']}"""
print(f'RUNNING {cmd}')
#os.system(cmd)

sample['arrbam']= os.path.join(out_qual,sample['name']+'_arr.bam')
cmd = f"""gatk ApplyBQSR \
            -R {ref_fasta} \
            -I {sample['markedbam']} \
            --bqsr-recal-file {sample['recaltable']} \
            -O {sample['arrbam']}"""
print(f'RUNNING {cmd}')
#os.system(cmd)

#variant calling
sample['vcf']= os.path.join(out_vcf,sample['name']+'.vcf')
sample['bamout']= os.path.join(out_vcf,sample['name']+'out.bam')
cmd = f"""gatk HaplotypeCaller \
            -R {ref_fasta} \
            -I {sample['arrbam']} \
            -O {sample['vcf']} \
            -bamout {sample['bamout']}"""
print(f'RUNNING {cmd}')
#os.system(cmd)

#CNN_filtering
#1D_CNN_filter
sample['1d_cnn_scored_vcf']= os.path.join(out_fil,sample['name']+'1d_cnn_scored.vcf')
cmd = f"""gatk CNNScoreVariants \
        -R {ref_fasta} \
        -V {sample['vcf']} \
        -O {sample['1d_cnn_scored_vcf']}"""
os.system(cmd)

sample['1d_cnn_filtered_vcf']= os.path.join(out_fil,sample['name']+'1d_cnn_filtered.vcf')
cmd = f"""gatk FilterVariantTranches \
        -V {sample['1d_cnn_scored_vcf']} \
        -O {sample['1d_cnn_filtered_vcf']} \
        --resource {resource} \
        --info-key CNN_1D \
        --snp-tranche 99.9 \
        --indel-tranche 95.0"""
#os.system(cmd)

#2D_CNN_filter
sample['2d_cnn_scored_vcf']= os.path.join(out_fil,sample['name']+'2d_cnn_scored.vcf')
cmd = f"""gatk CNNScoreVariants \
        -R {ref_fasta} \
        -I {sample['bamout']} \
        -V {sample['vcf']} \
        -O {sample['2d_cnn_scored_vcf']}
        --tensor-type read_tensor \
        --transfer-batch-size 8 \
        --inference-batch-size 8"""
#os.system(cmd)

sample['2d_cnn_filtered_vcf']= os.path.join(out_fil,sample['name']+'2d_cnn_filtered.vcf')
cmd = f"""gatk FilterVariantTranches \
        -V {sample['2d_cnn_scored_vcf']} \
        -O {sample['2d_cnn_filtered_vcf']} \
        --resource {resource} \
        --info-key CNN_2D \
        --snp-tranche 99.9 \
        --indel-tranche 95.0"""
#os.system(cmd)
