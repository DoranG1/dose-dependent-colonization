# Run filter2infer script on trimmed trimed reads.
rule filter2infer:
    params:
        outdir="workflow/out/{group}",
        trimdir="workflow/out/trimmed/{group}"
    output:
        seqtab="workflow/out/{group}/DADA2_output/seqtab_all.rds"
    conda:
        "Renv-minimal"
    script:
        "../scripts/filter2infer.R"

# Run makeTree script.
rule makeTree:
    input:
        "workflow/out/{group}/DADA2_output/seqtab_all.rds"
    params:
        outdir="workflow/out/{group}"
    output:
        tree="workflow/out/{group}/DADA2_output/dsvs_msa.tree"
    conda:
        "Renv-minimal"
    script:
        "../scripts/makeTree.R"

# Run summary_taxa script.
rule summary_taxa:
    input:
        "workflow/out/{group}/DADA2_output/seqtab_all.rds"
    params:
        outdir="workflow/out/{group}"
    output:
        "workflow/out/{group}/DADA2_output/ps_taxa.rds"
    conda:
        "Renv-minimal"
    script:
        "../scripts/summary_taxa.R"

# Run analyzeASVComposition script to generate a phyloseq object.
rule analyzeASVComposition:
    input:
        "workflow/out/{group}/DADA2_output/dsvs_msa.tree",
        "workflow/out/{group}/DADA2_output/ps_taxa.rds"
    params:
        outdir="workflow/out/{group}",
        metadata=config['samplesheet']
    output:
        ps="workflow/out/{group}/DADA2_output/ps_all.rds"
    conda:
        "Renv-minimal"
    script:
        "../scripts/analyzeASVComposition.R"
		
# Run convertPhyloseqToDataframe script to generate a dataframe.
rule convertPhyloseqToDataframe:
    input:
        "workflow/out/{group}/DADA2_output/ps_all.rds"
    params:
        outdir="workflow/out/{group}"
    output:
        ps="workflow/out/{group}/DADA2_output/ps_all.txt.gz"
    conda:
        "Renv-minimal"
    script:
        "../scripts/convertPhyloseqToDataframe.R"
