rule download_gtf:
    """ Download gtf from Ensembl """
    output:
        gtf="data/references/gtf/annotation.gtf"
    params:
        link=config["download"]["gtf"]
    shell:
        """
        mkdir -p data/references/gtf
        wget -O temp_gtf.gz {params.link}
        gunzip temp_gtf.gz
        mv temp_gtf {output.gtf}
        """

rule download_genome_fasta:
    """Download the reference genome from Ensembl."""
    output:
        genome="data/references/genome_fa/genome.fa"
    params:
        link=config["download"]["genome_fa"]
    shell:
        """
        mkdir -p data/references/genome_fa
        wget -O temp_genome.gz {params.link}
        gunzip temp_genome.gz
        mv temp_genome {output.genome}
        """

