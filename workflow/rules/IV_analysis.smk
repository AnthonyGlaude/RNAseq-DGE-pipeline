rule filter_discordant:
    """
    Sépare le BAM en deux : 
    - clean_bam (properly paired) pour featureCounts
    - discordant_bam (retirés) pour calcul du % et investigation future (splicing)
    """
    input:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam"
    output:
        clean_bam      = "results/STAR_filtered/{id}/Aligned.filtered.bam",
        discordant_bam = "results/qc_discordant/{id}/Aligned.discordant.bam"
    params:
        outdir_clean = "results/STAR_filtered/{id}",
        outdir_disc  = "results/qc_discordant/{id}"
    conda:
        "../envs/samtools.yml"
    threads:
        4
    log:
        "logs/{id}/filter_discordant.log"
    shell:
        """
        mkdir -p {params.outdir_clean} {params.outdir_disc} && \
        samtools view -@ {threads} -f 1 -F 2318 {input.bam} | cut -f1 | sort -u > {params.outdir_disc}/discordant_names.txt 2> {log} && \
        samtools view -@ {threads} -b -N ^{params.outdir_disc}/discordant_names.txt {input.bam} > {output.clean_bam} 2>> {log} && \
        samtools view -@ {threads} -b -N {params.outdir_disc}/discordant_names.txt {input.bam} > {output.discordant_bam} 2>> {log}
        """

CONDITIONS = ["Control_KO", "Control_WT", "HFpEF_KO", "HFpEF_WT"]

rule feature_counts:
    input:
        bams = lambda wc: expand(
            "results/STAR/{id}/Aligned.sortedByCoord.out.bam",
            id=[i for i in id_list if i.startswith(wc.condition)]
        ),
        gtf = rules.download_gtf.output.gtf
    output:
        counts  = "results/featurecounts/{condition}_counts.txt",
        summary = "results/featurecounts/{condition}_counts.txt.summary"
    params:
        outdir       = "results/featurecounts",
        strandedness = 2
    threads:
        8
    conda:
        "../envs/subread.yml"
    log:
        "logs/featurecounts/{condition}.log"
    shell:
        """
        mkdir -p {params.outdir} && \
        featureCounts \
        -a {input.gtf} \
        -o {output.counts} \
        -T {threads} \
        --tmpDir $SLURM_TMPDIR \
        -p --countReadPairs \
        -B -C \
        -M --fraction \
        -s {params.strandedness} \
        {input.bams} \
        &> {log}
        """

rule merge_counts:
    input:
        expand("results/featurecounts/{condition}_counts.txt", condition=CONDITIONS)
    output:
        merged="results/featurecounts/counts_merged.txt"
    shell:
        """
        module load r/4.3.1
        Rscript scripts/merge_counts.R {input} {output.merged}
        """

rule deseq2_analysis:
    """
    DESeq2 avec design condition * genotype :
    quatre comparaisons simples et un test d'interaction.
    """
    input:
        counts   = rules.merge_counts.output.merged,
        metadata = "data/references/sample_metadata.tsv"
    output:
        res_patho_WT = "results/deseq2/HFpEF_WT-Control_WT_DESeq2_gene.csv",
        res_KO_chow = "results/deseq2/Control_KO-Control_WT_DESeq2_gene.csv",
        res_interaction = "results/deseq2/condition_genotype_interaction_DESeq2_gene.csv",
        res_patho_KO = "results/deseq2/HFpEF_KO-Control_KO_DESeq2_gene.csv",
        res_KO_patho = "results/deseq2/HFpEF_KO-HFpEF_WT_DESeq2_gene.csv",
        dds = "results/deseq2/dds.rds"
    params:
        outdir = "results/deseq2",
        filter_count_threshold = config["dge"]["filter_count_threshold"]
    conda:
        "../envs/DESeq2.yml"
    log:
        "logs/deseq2/analysis.log"
    script:
        "../scripts/deseq2_analysis.R"



rule deseq2_stats:
    """
    Analyses post-DESeq2 pour chaque contraste :
    tableau annoté, volcano plots et enrichissements GO.
    """
    input:
        deseq2 = "results/deseq2/{comp}_DESeq2_gene.csv"
    output:
        stat = "results/deseq2_stats/{comp}_stats.csv",
        volcano_total = "results/volcano/{comp}_volcano_total.png",
        volcano_zoom = "results/volcano/{comp}_volcano_zoom.png",

        deg_total = "results/deg_lists/{comp}_all_DEG.tsv",
        deg_up = "results/deg_lists/{comp}_upregulated.tsv",
        deg_down = "results/deg_lists/{comp}_downregulated.tsv",

        go_total = "results/go/{comp}_GO_all.png",
        go_enrich_up = "results/go/{comp}_GO_up.png",
        go_enrich_down = "results/go/{comp}_GO_down.png"
    log:
        "logs/deseq2/stats_{comp}.log"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/deseq2_stats.R"


rule PCA_MA_Scatter_plot_deseq2:
    input:
        dds = rule.deseq2_analysis.output.dds
    output:
        pca = "results/figure/total_pca.png"
        scatter = "results/figure/{comp}_scatter.png"
        ma = "results/figure/{comp}_ma.png"
    log:
        "logs/deseq2/figures_{comp}.log"
    conda:
        "../envs/DESeq2.yml"
    script:
        "../scripts/figures.R"


rule discordant_reads:
    """ Calcule du % de reads discordants (flag -F 1294) par échantillon """
    input:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam"
    output:
        stats = "results/qc_discordant/{id}_discordant.txt"
    conda:
        "../envs/samtools.yml"
    log:
        "logs/{id}/discordant.log"
    shell:
        """
        total=$(samtools view -c -f 1 {input.bam})
        discordant=$(samtools view -c -F 1294 {input.bam})
        pct=$(echo "scale=4; $discordant / $total * 100" | bc)
        echo -e "{wildcards.id}\\t$total\\t$discordant\\t$pct" > {output.stats}
        """

rule aggregate_discordant:
    """ Rassemble les stats de discordance de tous les échantillons en un seul tableau """
    input:
        stats = expand("results/qc_discordant/{id}_discordant.txt", id=id_list)
    output:
        table = "results/qc_discordant/all_samples_discordant.tsv"
    shell:
        """
        echo -e "sample_id\\ttotal_reads\\tdiscordant_reads\\tpct_discordant" > {output.table}
        cat {input.stats} >> {output.table}
        """

rule discordant_qc_plot:
    """
    Boxplot du % de reads discordants par groupe (chow_WT, chow_KO, patho_WT, patho_KO)
    + test Kruskal-Wallis pour vérifier l'absence de différence significative
    entre groupes (contrôle de biais technique).
    """
    input:
        table    = rules.aggregate_discordant.output.table,
        metadata = "/home/glaudea/scratch/glaudea/test_souris_rnaseq/RNA_seq-analysis/workflow/data/references/sample_metadata.tsv"
    output:
        plot = "results/qc_discordant/discordant_boxplot.svg",
        test = "results/qc_discordant/kruskal_test_result.txt"
    conda:
        "../envs/DESeq2.yml"
    log:
        "logs/qc_discordant/plot.log"
    script:
        "../scripts/discordant_qc_plot.R"


rule feature_counts_test_correction:
    input:
        bams = [
            "results/STAR/Control_KO_25/Aligned.sortedByCoord.out.bam",
        ],
        gtf = rules.download_gtf.output.gtf
    output:
        counts  = "results/featurecounts_test/Control_KO_counts.txt",
        summary = "results/featurecounts_test/Control_KO_counts.txt.summary"
    params:
        outdir = "results/featurecounts_test"
    threads:
        8
    conda:
        "../envs/subread.yml"
    log:
        "logs/featurecounts/Control_KO_test_correction.log"
    shell:
        """
        mkdir -p {params.outdir} && \
        featureCounts \
        -a {input.gtf} \
        -o {output.counts} \
        -T {threads} \
        --tmpDir $SLURM_TMPDIR \
        -p --countReadPairs \
        -B -C -P \
        -O -M --fraction \
        -s 2 \
        {input.bams} \
        &> {log}
        """
