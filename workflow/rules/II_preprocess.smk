# RNA-Seq quality control followed by alignment with STAR 
rule fastqc:
    """Assess the FASTQ quality using FastQC BEFORE TRIMMING"""
    input:
        fq1 = os.path.join(config["path"]["fastq"], "{id}_R1.fastq.gz"),
        fq2 = os.path.join(config["path"]["fastq"], "{id}_R2.fastq.gz")
    output:
        html1 = "data/qc/{id}/{id}_R1_fastqc.html",
        html2 = "data/qc/{id}/{id}_R2_fastqc.html",
    params:
        out_dir = "data/qc/{id}"
    log:
        "logs/{id}/fastqc.log"
    threads: 8
    conda:
        "../envs/fastqc.yml" 
    shell:
        "mkdir -p {params.out_dir} && "
        "fastqc --outdir {params.out_dir} --format fastq --threads {threads} {input.fq1} {input.fq2} &> {log} "

rule trim_reads:
    input:
        fq1 = os.path.join(config["path"]["fastq"], "{id}_R1.fastq.gz"),
        fq2 = os.path.join(config["path"]["fastq"], "{id}_R2.fastq.gz")
    output:
        trim1     = "data/trim_galore/{id}/{id}_R1_val_1.fq.gz",
        trim2     = "data/trim_galore/{id}/{id}_R2_val_2.fq.gz",
        unpaired1 = "data/trim_galore/{id}/{id}_R1_unpaired_1.fq.gz",
        unpaired2 = "data/trim_galore/{id}/{id}_R2_unpaired_2.fq.gz"
    params:
        out_dir = "data/trim_galore/{id}"
    threads:
        6
    conda:
        "../envs/trim_galore.yml"
    log:
        "logs/{id}/trim.log"
    shell:
        """
        mkdir -p {params.out_dir} &&\
        trim_galore --paired \
        --retain_unpaired \
        --cores {threads} \
        --gzip \
        --output_dir {params.out_dir} \
        {input.fq1} {input.fq2} \
        &> {log}
        """

# FASTP trimming
#rule trim_reads:
#    input:
#        fq1 = os.path.join(config["path"]["fastq"], "{id}_R1.fastq.gz"),
#        fq2 = os.path.join(config["path"]["fastq"], "{id}_R2.fastq.gz")
#    output:
#        trim1     = "data/fastp/{id}/{id}_R1_trimmed.fastq.gz",
#        trim2     = "data/fastp/{id}/{id}_R2_trimmed.fastq.gz",
#        unpaired1 = "data/fastp/{id}/{id}_R1_unpaired.fastq.gz",
#        unpaired2 = "data/fastp/{id}/{id}_R2_unpaired.fastq.gz",
#        html      = "data/fastp/{id}/{id}_fastp.html",
#        json      = "data/fastp/{id}/{id}_fastp.json"
#    params:
#        out_dir = "data/fastp/{id}",
#        options = "--qualified_quality_phred 20 "
#                  "--length_required 20 "
#                  "--cut_window_size 4 "
#                  "--cut_mean_quality 20 "
#                  "--cut_front "
#                  "--cut_tail "
#                  "--low_complexity_filter "
#                  "--complexity_threshold 30"
#    threads:
#        6
#    conda:
#        "../envs/fastp.yml"
#    log:
#        "logs/{id}/trim.log"
#    shell:
#        """
#        mkdir -p {params.out_dir} &&\
#        fastp \
#          -i {input.fq1} -I {input.fq2} \
#          -o {output.trim1} -O {output.trim2} \
#          --unpaired1 {output.unpaired1} \
#          --unpaired2 {output.unpaired2} \
#          --thread {threads} \
#          -h {output.html} \
#          -j {output.json} \
#          {params.options} \
#        &> {log}
#        """


rule qc_fastq:
    """ Assess the FASTQ quality using FastQC AFTER TRIMMING"""
    input:
        trimm_fq1 = rules.trim_reads.output.trim1,
        trimm_fq2 = rules.trim_reads.output.trim2,
        unpaired1 = rules.trim_reads.output.unpaired1,
        unpaired2 = rules.trim_reads.output.unpaired2
    output:
        qc_trimm_fq1_out = "data/qc_after_trim/{id}/{id}_R1_val_1_fastqc.html",
        qc_trimm_fq2_out = "data/qc_after_trim/{id}/{id}_R2_val_2_fastqc.html"
    params:
        out_dir = "data/qc_after_trim/{id}"
    log:
        "logs/{id}/FASTQC2.log"
    threads:
        8
    conda:
        "../envs/fastqc.yml"
    shell:
        """
        mkdir -p {params.out_dir} &&
        fastqc \
            --outdir {params.out_dir} \
            --format fastq \
            --threads {threads} \
            {input.trimm_fq1} {input.trimm_fq2} {input.unpaired1} {input.unpaired2}\
            &> {log}
        """
