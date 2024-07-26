rule copy_reference:
    input:
        fasta = get_reference_fasta
    output:
        "resources/{ref}/{ref}.fasta"
    shell:
        """
        cp {input.fasta} {output}
        """

checkpoint genome_faidx:
    input:
        path = rules.copy_reference.output
    output:
        "resources/{ref}/{ref}.fasta.fai"
    cache: True
    conda:
        "ngs"
    shell:
        """
        samtools faidx {input.path}
        """
        

