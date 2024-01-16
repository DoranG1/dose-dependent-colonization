# Run filter2infer script on trimmed trimed reads.
rule filter2infer:
    params:
        outdir="workflow/out/{group}",
        trimdir="workflow/out/trimmed/{group}"
    output:
        seqtab="workflow/out/{group}/DADA2_output/seqtab_all.rds"
    conda:
        "../envs/snakemakeR4-0.yml"
    script:
        "../scripts/filter2infer.R"

# Run makeTree script.
rule makeTree:
    params:
        outdir="workflow/out/{group}"
    output:
        tree="workflow/out/{group}/DADA2_output/dsvs_msa.tree"
    conda:
        "../envs/snakemakeR4-0.yml"
    script:
        "../scripts/makeTree.R"

# Run summary_taxa script.
rule summary_taxa:
    params:
        outdir="workflow/out/{group}"
    output:
        summary="workflow/out/{group}/DADA2_output/L6_summary.txt"
    conda:
        "../envs/snakemakeR4-0.yml"
    script:
        "../scripts/summary_taxa.R"

# Run analyzeASVComposition script.
rule analyzeASVComposition:
    params:
        outdir="workflow/out/{group}"
    output:
        ps="workflow/out/{group}/DADA2_output/ps_all.rds"
    conda:
        "../envs/snakemakeR4-0.yml"
    script:
        "../scripts/analyzeASVComposition.R"
