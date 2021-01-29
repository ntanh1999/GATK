#single sample
#from fastq

import os

out_dir = "/media/vinbdi/data/tienanh/gatk"
os.mkdir(out_dir)
out_aln = os.path.join(out_dir,'aln')
os.mkdir(out_aln)
out_qual = os.path.join(out_dir,'qualification')
os.mkdir(out_qual)
tmp = os.path.join(out_dir,'tmp')
os.mkdir(tmp)
out_vcf = os.path.join(out_dir,'vcf')
os.mkdir(out_vcf)

ref_fasta = '/home/vinbdi/Desktop/ref38/hg38.fasta'
dbsnp = '/home/vinbdi/Desktop/ref38/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf'

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
    group['bam']= os.path.join(out_aln,sample['name']+'_'+group['groupname']+'_mapped.bam')
    cmd = f"""bwa mem -M \
                -t 20 \
                {ref_fasta} \
                {group['read1']} \
                {group['read2']} \
                | samtools view -Shb -o {group['bam']}"""
    print(f'RUNNING {cmd}')
    os.system(cmd)
    
    #add readgroup metadata
    group['mbam']= os.path.join(out_aln,sample['name']+'_'+group['groupname']+'_2_mapped.bam')
    cmd = f"""gatk AddOrReplaceReadGroups \
            -I {group['bam']} \
            -O {group['mbam']} \
            --RGID {group['groupname']} \
            --RGSM {sample['name']} \
            --RGLB lib1 \
            --RGPU unit1 \
            --RGPL ILLUMINA"""
    print(f'RUNNING {cmd}')
    os.system(cmd)
            
    #sortbam
    group['sortedbam']= os.path.join(out_aln,sample['name']+'_'+group['groupname']+'_sorted.bam')
    cmd = f"""gatk SortSam \
            -I {group['mbam']} \
            -O {group['sortedbam']} \
            -SORT_ORDER coordinate \
            --TMP_DIR {tmp}"""
    print(f'RUNNING {cmd}')
    os.system(cmd)

#mark duplicate
group['markedbam']= os.path.join(out_qual,sample['name']+'_marked.bam')
group['metrics']= os.path.join(out_qual,sample['name']+'_metrics.txt')

cmd = f"""gatk MarkDuplicates \
            -I {sample['groups'][0]['sortedbam']} \
            -I {sample['groups'][1]['sortedbam']} \
            -I {sample['groups'][2]['sortedbam']} \
            -I {sample['groups'][3]['sortedbam']} \
            -O {group['markedbam']} \
            -M {group['metrics']} \
            --TMP_DIR {tmp}"""
print(f'RUNNING {cmd}')
os.system(cmd)

#base calibration
group['recaltable']= os.path.join(out_qual,sample['name']+'_recal.table')
cmd = f"""gatk BaseRecalibrator \
        -I {group['markedbam']} \
        -R {ref_fasta} \
        --known-sites {dbsnp} \
        -O {group['recaltable']}"""
print(f'RUNNING {cmd}')
os.system(cmd)

group['arrbam']= os.path.join(out_qual,sample['name']+'_arr.bam')
cmd = f"""gatk ApplyBQSR \
            -R {ref_fasta} \
            -I {group['markedbam']} \
            --bqsr-recal-file {group['recaltable']} \
            -O {group['arrbam']}"""
print(f'RUNNING {cmd}')
#os.system(cmd)

#variant calling
group['vcf']= os.path.join(out_vcf,sample['name']+'.vcf')
cmd = f"""gatk HaplotypeCaller \
            -R {ref_fasta} \
            -I {group['arrbam']} \
            -O {group['vcf']}"""
print(f'RUNNING {cmd}')
os.system(cmd)