rule merge_variants:
    input:
        vcfs = get_sample_first_vcfs,
        seq_ref_list = "resources/{ref}/{ref}.fasta.fai",
    params:
        mem = "-Xmx300g",
        
    resources:
         mem_mb=40000
    output:
        f'{base_dir}/results/mapping/{{ref}}/vcf/merged_variants.vcf'
    log:
        'log/mapping/{ref}/vcf/merge_variants.log'
    conda:
        'ngs'
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            MergeVariants -s {input.seq_ref_list} -o {output} {input.vcfs} 2> {log}
        """

rule merge_bams:
    input: 
        get_bams_to_merge_with_samples
    output:
        f'{base_dir}/results/mapping/{{ref}}/duplicates/{{sample}}_merged.bam'
    conda:
        "ngs"
    shell:
        """
        samtools merge {output} {input}
        """

rule single_sample_variant_detector_two:
    input:
        bam = get_bam_variant_detection_two,
        known_variants = f'{base_dir}/results/mapping/{{ref}}/vcf/merged_variants.vcf',
        ref = "resources/{ref}/{ref}.fasta"
    output:
        vcf = f"{base_dir}/results/mapping/{{ref}}/vcf/individual/{{sample}}_NGSEP.vcf.gz"
    params:
        params = " ".join(["-{param} {value}".format(param=i_param, value = config['NGSEP']['SingleSampleVariantsDetector'][i_param]) for i_param in config['NGSEP']['SingleSampleVariantsDetector'].keys()]),
        mem = "-Xmx3g",
    resources:
         mem_mb=3000
    log:
        "log/mapping/{ref}/vcf/individual/{sample}_NGSEP.log"
    conda:
        "ngs"
    shell:
        """
        java {params.mem} -jar {config[NGSEP][path]} \
            SingleSampleVariantsDetector {params.params} -knownVariants {input.known_variants}\
            -sampleId {wildcards.sample} \
            -r {input.ref} \
            -i {input.bam} \
            -o {output.vcf} 2> {log} &&
            cat {output.vcf}.vcf| bgzip > {output.vcf} &&
            rm {output.vcf}.vcf
        """

rule index_vcfs:
    input:
        vcf = f"{base_dir}/results/mapping/{{ref}}/vcf/individual/{{sample}}_NGSEP.vcf.gz"
    output:
        index = f"{base_dir}/results/mapping/{{ref}}/vcf/individual/{{sample}}_NGSEP.vcf.gz.tbi"
    conda:
        "ngs"
    shell:
        """
        if [[ ! -f "{input.vcf}.tbi" ]]
        then
            tabix -p vcf {input.vcf}
        fi
        """

rule split_vcfs_by_chr:
    input:
        vcf = f"{base_dir}/results/mapping/{{ref}}/vcf/individual/{{sample}}_NGSEP.vcf.gz"
    output:
        region_vcf = temp(f"{base_dir}/results/mapping/{{ref}}/vcf/concat/{{chrom}}/{{sample}}_{{chrom}}_NGSEP.vcf.gz"),
        region_vcf_idx = temp(f"{base_dir}/results/mapping/{{ref}}/vcf/concat/{{chrom}}/{{sample}}_{{chrom}}_NGSEP.vcf.gz.tbi")
    conda:
        "ngs"
    shell:
        """
        bcftools view -r {wildcards.chrom} {input.vcf}| bgzip > {output.region_vcf} && \
        tabix -p vcf {output.region_vcf}
        """

rule concat_by_chr:
    input:
        unpack(get_by_chr_vcfs)
    output:
        vcf = temp(f"{base_dir}/results/mapping/{{ref}}/vcf/concat/{{chrom}}/{{chrom}}_NGSEP.vcf.gz"),
        index = temp(f"{base_dir}/results/mapping/{{ref}}/vcf/concat/{{chrom}}/{{chrom}}_NGSEP.vcf.gz.tbi"),
    conda:
        "ngs"
    shell:
        """
        bcftools concat {input.vcfs}| bgzip > {output.vcf} && \
        tabix -p vcf {output.vcf}
        """

rule concat_all:
    input:
        unpack(get_chr_vcfs)
    output:
        vcf = f"{base_dir}/results/mapping/{{ref}}/vcf/concat/merged_NGSEP.vcf.gz",
        index = f"{base_dir}/results/mapping/{{ref}}/vcf/concat/merged_NGSEP.vcf.gz.tbi",
    conda:
        "ngs"
    shell:
        """
        bcftools concat {input.vcfs}| bgzip > {output.vcf} && \
        tabix -p vcf {output.vcf}
        """
