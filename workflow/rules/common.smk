import pandas as pd
import os
import glob
from yaml import safe_load
import random
import string
import numpy as np


## Configuration setup        
configfile: "config/config.yaml"

base_dir = config['base_dir']
## Read population sample file
samples = pd.read_csv(config['samples_file'], sep='\t')
samples['sample_name'] = samples['line_id']

duplicates = samples[samples.duplicated(subset = 'line_id',keep=False)].copy()
duplicates['new_name'] = duplicates.line_id + '-' + duplicates.plate.astype(str)




duplicate_dict = pd.Series(duplicates.line_id.values,index=duplicates.new_name).to_dict()
dups = [x[0] for x in duplicate_dict.items()]

#print(duplicate_dict)
#print(dups)


samples.loc[duplicates.index, 'line_id'] = duplicates['new_name'].tolist()
samples.set_index('line_id', drop=False, inplace=True)

references = pd.read_table(config['references'], sep = '\t')
references.set_index("ref_name", inplace = True)

def get_reference_fasta(wildcards):
    print(references.loc[wildcards.ref, "ref_path"])
    return(references.loc[wildcards.ref, "ref_path"])

def get_sample_first_vcfs(wildcards):
    paths = list()
    for n, row in samples.iterrows():
        sample_name = row['sample_name']
        vcf_path = row['path'] + "_bowtie2_NGSEP.vcf.gz" 
        paths.append(vcf_path)
    return paths


def get_bams_to_merge_with_samples(wildcards):
    target_samples = samples[samples['sample_name'] == wildcards.sample]
    paths = list()
    for i,row in target_samples.iterrows():
        paths.append(row['path'] + "_bowtie2_sorted.bam")
    return paths

def get_bam_variant_detection_two(wildcards): 
    try:
        if wildcards.sample in list(duplicate_dict.keys()):
            sample_name = duplicate_dict[wildcards.sample]
        else:
            sample_name = wildcards.sample
        sample_data = samples.loc[wildcards.sample]
        bam_path = sample_data['path'] + '_bowtie2_sorted.bam'
        
    except Exception as e:
        print(e)
        bam_path = "{base_dir}/results/mapping/{ref}/duplicates/{sample}_merged.bam".format(base_dir = base_dir, ref = wildcards.ref, sample = sample_name)
    return bam_path 

def get_by_chr_vcfs(wildcards):    
    chroms = ['Chr01', 'Chr02', 'Chr03', 'Chr04', 'Chr05', 'Chr06','Chr07','Chr08','Chr09','Chr10', 'Chr11']
    samples_unique = samples['sample_name'].unique()
    vcfs_list = list()
    index_list = list()
    for chrom in chroms:
        for sample in samples_unique:
        vcfs = ["{base_dir}/results/mapping/{ref}/vcf/concat/{chrom}/{sample}_{chrom}_NGSEP.vcf.gz".format(
            base_dir = base_dir,ref = wildcards.ref, sample = x, chrom = chrom) for x in samples_unique]
        index = ["{base_dir}/results/mapping/{ref}/vcf/concat/{chrom}/{sample}_{chrom}_NGSEP.vcf.gz.tbi".format(
            base_dir = base_dir,ref = wildcards.ref, sample = x, chrom = chrom) for x in samples_unique]
        vcfs_list.extend(vcfs)
        index_list.extend(index)
    out = {'vcfs': vcfs_list, 'indices': index_list}
    return out

def get_chr_vcfs(wildcards):    
    chroms = ['Chr01', 'Chr02', 'Chr03', 'Chr04', 'Chr05', 'Chr06','Chr07','Chr08','Chr09','Chr10', 'Chr11']

    vcfs = ["{base_dir}/results/mapping/{ref}/vcf/concat/{chrom}/{chrom}_NGSEP.vcf.gz".format(
        base_dir = base_dir,ref = wildcards.ref, chrom = x) for x in chroms]
    index = ["{base_dir}/results/mapping/{ref}/vcf/concat/{chrom}/{chrom}_NGSEP.vcf.gz.tbi".format(
        base_dir = base_dir,ref = wildcards.ref, chrom = x) for x in chroms]

    out = {'vcfs': vcfs, 'indices': index}
    return out



def get_population_vcfs(wildcards):    
    
    samples_unique = samples['sample_name'].unique()
    vcfs = ["{base_dir}/results/mapping/{ref}/vcf/individual/{sample}_NGSEP.vcf.gz".format(
        base_dir = base_dir,ref = wildcards.ref, sample = x) for x in samples_unique]
    index = ["{base_dir}/results/mapping/{ref}/vcf/individual/{sample}_NGSEP.vcf.gz.tbi".format(
        base_dir = base_dir, ref = wildcards.ref, sample = x) for x in samples_unique]
    out = {'vcfs': vcfs, 'index': index}
    return out

wildcard_constraints:
    ref = "|".join(references.index),

