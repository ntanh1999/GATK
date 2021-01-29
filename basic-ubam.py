#single sample
#convert to ubam and merge

import os

out_dir = "/home/ted/ubuntu/ADR/out"
#os.mkdir(out_dir)

ref_fasta = 'path_ref_fasta'
dbsnp = 'path_dbsnp'

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
    #convert fastq to ubam
    group['ubam']=os.path.join(out_aln,sample['name']+'_'+group['groupname']+'_unmapped.bam')
    cmd =  f"""gatk FastqToSam \
                -F1 {group['read1']} \
                -F2 {group['read2']} \
                -O {group['ubam']} \
                -SM {sample['name']} \
                -RG {group['groupname']} \
                -LB lib1 \
                -PL ILLUMINA"""
    #os.system(cmd)

    #map reads to reference
    group['mbam']= os.path.join(out_aln,sample['name']+'_'+group['groupname']+'_mapped.bam')
    cmd = f"""bwa mem -M \
                -t 20 \
                -R "@RG\tID:{group['groupname']}\tSM:{sample['name']}\tLB:lib1\tPL:ILLUMINA" \
                {ref_fasta} \
                {group['read1']} \
                {group['read2']} \
                | samtools view -Shb -o {group['mbam']}"""
    #os.system(cmd)

    #merge bam alignment
    group['bothbam']= os.path.join(out_aln,sample['name']+'_'+group['groupname']+'_both.bam')
    cmd = f"""gatk MergeBamAlignment \
                -ALIGNED {group['mbam']} \
                -UNMAPPED {group['ubam']} \
                -O {group['bothbam']} \
                -R {ref_fasta}"""
    #os.system(cmd)

#merge, mark duplicate, sort

out_qual = os.path.join(out_dir,'qualification')
#os.mkdir(out_qual)

tmp = os.path.join(out_dir,'tmp')
#os.mkdir(tmp)

group['markedbam']= os.path.join(out_qual,sample['name']+'_marked.bam')
group['metrics']= os.path.join(out_qual,sample['name']+'_metrics.txt')

cmd = f"""gatk MarkDuplicates \
            -I {sample['groups'][0]['bothbam']} \
            -I {sample['groups'][1]['bothbam']} \
            -I {sample['groups'][2]['bothbam']} \
            -I {sample['groups'][3]['bothbam']} \
            -O {group['markedbam']} \
            -M {group['metrics']} \
            --TMP_DIR {tmp}"""
#os.system(cmd)

group['sortedbam']= os.path.join(out_qual,sample['name']+'_sorted.bam')
cmd = f"""gatk SortSam \
        -I {group['markedbam']} \
        -O {group['sortedbam']} \
        -SORT_ORDER coordinate \
        --TMP_DIR {tmp}"""
#os.system(cmd)

group['index']= os.path.join(out_qual,sample['name']+'_')
cmd = f'samtools index {group['sortedbam']}'
#os.system(cmd)

#base calibration
group['recaltable']= os.path.join(out_qual,sample['name']+'_recal.table')
cmd = f"""gatk BaseRecalibrator \
        -I {group['sortedbam']} \
        -R {ref_fasta} \
        --known-sites {dbsnp} \
        -O {group['recaltable']}"""
#os.system(cmd)

group['arrbam']= os.path.join(out_qual,sample['name']+'_arr.bam')
cmd = f"""gatk ApplyBQSR \
            -R {ref_fasta} \
            -I {group['sortedbam']} \
            --bqsr-recal-file {group['recaltable']} \
            -O {group['arrbam']}"""
#os.system(cmd)

#variant calling
out_vcf = os.path.join(out_dir,'vcf')
#os.mkdir(out_vcf)

group['vcf']= os.path.join(out_qual,sample['name']+'.vcf')
cmd = f"""gatk HaplotypeCaller \
            -R {ref_fasta} \
            -I {group['arrbam']} \
            -O {group['vcf']}"""
#os.system(cmd)


print(sample)
