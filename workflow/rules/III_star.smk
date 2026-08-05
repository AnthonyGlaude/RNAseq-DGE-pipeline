rule star_index:
    """ Generates the genome index for STAR """
    input:
        fasta = rules.download_genome_fasta.output.genome,
        gtf = rules.download_gtf.output.gtf # reference
    output:
        chrNameLength = "data/references/star_index/chrNameLength.txt"
    params:
        dir = config['path']['star_index']  
    log:
        "logs/index.log"
    conda:
        "../envs/star.yml"
    threads:
        32
    shell:
        """
        mkdir -p {params.dir} && \
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params.dir} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 79 \
        &> {log}
        """

rule star_alignreads:
    """ Generates a bam file using STAR
        parameters based on Sebastien Leblanc and Marie A. Brunet, 2020 methodology,
        adapted for 75nt reads
    """
    input:
        idx = rules.star_index.output,
        fq1 = rules.trim_reads.output.trim1,
        fq2 = rules.trim_reads.output.trim2
    output:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam",
        bam_logs = "results/STAR/{id}/Log.final.out"
    params:
        index = config['path']['star_index'],
        output_dir = "results/STAR/{id}/"
    log:
        "logs/{id}/alignreads.log"
    threads:
        32
    conda:
        "../envs/star.yml"
    shell:
        """
        mkdir -p {params.output_dir} && \
        STAR --runMode alignReads \
            --genomeDir {params.index} \
            --readFilesIn {input.fq1} {input.fq2} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outReadsUnmapped Fastx \
            --outFilterType BySJout \
            --outStd Log \
            --outSAMunmapped None \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.output_dir} \
            --outFilterMismatchNmax 5 \
            --outFilterMultimapNmax 10 \
            --outSAMprimaryFlag AllBestScore \
            --limitBAMsortRAM 90000000000 \
            --outTmpDir /tmp/{wildcards.id} \
            &> {log}
        """
