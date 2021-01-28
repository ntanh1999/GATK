#multi sample
#from fastq

import os
import pandas as pd

out_dir = "/home/ted/ubuntu/ADR/out"
#os.mkdir(out_dir)
out_aln = os.path.join(out_dir,'aln')
#os.mkdir(out_aln)
out_qual = os.path.join(out_dir,'qualification')
#os.mkdir(out_qual)
tmp = os.path.join(out_dir,'tmp')
#os.mkdir(tmp)
out_vcf = os.path.join(out_dir,'vcf')
#os.mkdir(out_vcf)

ref_fasta = 'path_ref_fasta'
dbsnp = 'path_dbsnp'

input = pd.read_csv("input.tsv",sep='\t', dtype='str')
input.fillna('', inplace=True)

collection={}
collection['samples'] = {}
for i,row in input.iterrows():
    if row['name'] not in collection['samples'].keys():
        collection['samples'][row['name']]={}
    collection['samples'][row['name']][row['group']]={} 
    collection['samples'][row['name']][row['group']]['read1']=row['read1']
    collection['samples'][row['name']][row['group']]['read2']=row['read2']

list_sample = []
for sample_name in collection['samples']:
    list_sample.append(sample_name)
    sample=collection['samples'][sample_name]
    list_group = []
    for group_name in sample:
        group = sample[group_name]
        list_group.append(group_name)

        #map reads to reference
        group['mappedbam']= os.path.join(out_aln,sample_name+'_'+group_name+'_mapped.bam')
        cmd = f"""bwa mem -M \
                    -t 20 \
                    -R @RG\tID:${group_name}\tSM:${sample_name}\tLB:lib1\tPL:ILLUMINA \
                    {ref_fasta} \
                    {group['read1']} \
                    {group['read2']} \
                    | samtools view -Shb -o {group['mappedbam']}"""
        #os.system(cmd)

    #mark duplicate
    sample['markedbam']= os.path.join(out_qual,sample_name+'_marked.bam')
    sample['metrics']= os.path.join(out_qual,sample_name+'_metrics.txt')

    cmd = f"""gatk MarkDuplicates \
                -O {sample['markedbam']} \
                -M {sample['metrics']} \
                --TMP_DIR {tmp}"""
    for gn in list_group:
        cmd += f"-I {sample[gn]['mappedbam']}"

    #os.system(cmd)

    #base calibration
    sample['recaltable']= os.path.join(out_qual,sample_name+'_recal.table')
    cmd = f"""gatk BaseRecalibrator \
            -I {sample['markedbam']} \
            -R {ref_fasta} \
            --known-sites {dbsnp} \
            -O {sample['recaltable']}"""
    #os.system(cmd)
    sample['arrbam']= os.path.join(out_qual,sample_name+'_arr.bam')
    cmd = f"""gatk ApplyBQSR \
                -R {ref_fasta} \
                -I {sample['markedbam']} \
                --bqsr-recal-file {sample['recaltable']} \
                -O {sample['arrbam']}"""
    #os.system(cmd)

    #variant calling
    sample['gvcf']= os.path.join(out_vcf,sample_name+'_g.vcf.gz')
    cmd = f"""gatk HaplotypeCaller \
                -R {ref_fasta} \
                -I {sample['arrbam']} \
                -O {sample['gvcf']} \
                -ERC GVCF"""
    #os.system(cmd)

#consolidating GVCFs
out_joint = os.path.join(out_dir,'joint')
#os.mkdir(out_joint)
out_gdb = os.path.join(out_joint,'genomicsdb')
#os.mkdir(out_gdb)
collection['genomicsdb'] = out_gdb
cmd = f"""gatk GenomicsDBImport \
        --genomicsdb-workspace-path {collection['genomicsdb']}"""
for i in list_sample:
    cmd += f"-V {collection['samples'][i]['gvcf']}"
#os.system(cmd)

collection['jointvcf'] = os.path.join(out_joint,'joint.vcf.gz')
cmd =f"""gatk GenotypeGVCFs
        -R {ref_fasta} \
        -V {collection['genomicsdb']} \
        -O {collection['jointvcf']}"""
#os.system(cmd) 

print(collection)