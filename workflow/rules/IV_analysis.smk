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


#######################Reanalysis with wee1-as integrated in annotations ######################

rule feature_counts_wee1as:
    input:
        bams = lambda wc: expand(
            "results/STAR/{id}/Aligned.sortedByCoord.out.bam",
            id=[i for i in id_list if i.startswith(wc.condition)]
        ),
        gtf = "/home/glaudea/scratch/glaudea/test_souris_rnaseq/RNA_seq-analysis/workflow/data/references/gtf/annotation_wee1_as.gtf"

    output:
        counts  = "results_avec_wee1as/featurecounts/{condition}_counts.txt",
        summary = "results_avec_wee1as/featurecounts/{condition}_counts.txt.summary"

    params:
        outdir = "results_avec_wee1as/featurecounts",
        strandedness = 2

    threads:
        8

    conda:
        "../envs/subread.yml"

    log:
        "logs/featurecounts_wee1as/{condition}.log"

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

rule merge_counts_wee1as:
    input:
        expand(
            "results_avec_wee1as/featurecounts/{condition}_counts.txt",
            condition=CONDITIONS
        )
    output:
        merged = "results_avec_wee1as/featurecounts/counts_merged.txt"
    shell:
        """
        module load r/4.3.1
        Rscript scripts/merge_counts.R {input} {output.merged}
        """

rule deseq2_analysis_wee1as:
    input:
        counts = rules.merge_counts_wee1as.output.merged,
        metadata = "data/references/sample_metadata.tsv"

    output:
        res_patho_WT = "results_avec_wee1as/deseq2/HFpEF_WT-Control_WT_DESeq2_gene.csv",
        res_KO_chow = "results_avec_wee1as/deseq2/Control_KO-Control_WT_DESeq2_gene.csv",
        res_interaction = "results_avec_wee1as/deseq2/condition_genotype_interaction_DESeq2_gene.csv",
        res_patho_KO = "results_avec_wee1as/deseq2/HFpEF_KO-Control_KO_DESeq2_gene.csv",
        res_KO_patho = "results_avec_wee1as/deseq2/HFpEF_KO-HFpEF_WT_DESeq2_gene.csv",
        dds = "results_avec_wee1as/deseq2/dds.rds"

    params:
        outdir = "results_avec_wee1as/deseq2",
        filter_count_threshold = config["dge"]["filter_count_threshold"]

    conda:
        "../envs/DESeq2.yml"

    log:
        "logs/deseq2_wee1as/analysis.log"

    script:
        "../scripts/deseq2_analysis.R"

rule deseq2_stats_wee1as:
    input:
        deseq2 = "results_avec_wee1as/deseq2/{comp}_DESeq2_gene.csv"

    output:
        stat = "results_avec_wee1as/deseq2_stats/{comp}_stats.csv",

        volcano_total = "results_avec_wee1as/volcano/{comp}_volcano_total.png",
        volcano_zoom = "results_avec_wee1as/volcano/{comp}_volcano_zoom.png",

        deg_total = "results_avec_wee1as/deg_lists/{comp}_all_DEG.tsv",
        deg_up = "results_avec_wee1as/deg_lists/{comp}_upregulated.tsv",
        deg_down = "results_avec_wee1as/deg_lists/{comp}_downregulated.tsv",

        go_total = "results_avec_wee1as/go/{comp}_GO_all.png",
        go_enrich_up = "results_avec_wee1as/go/{comp}_GO_up.png",
        go_enrich_down = "results_avec_wee1as/go/{comp}_GO_down.png"

    log:
        "logs/deseq2_wee1as/stats_{comp}.log"

    conda:
        "../envs/DESeq2.yml"

    script:
        "../scripts/deseq2_stats.R"
