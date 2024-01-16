# Use umi-tools to annotate in-line barcodes
# based on the first 7 bases of each sequencing read.
rule demultiplexBarcodes:
    input:
        r1= lambda wildcards: dffastq.query('file == @wildcards.sample')['read1'].tolist()[0],
        r2= lambda wildcards: dffastq.query('file == @wildcards.sample')['read2'].tolist()[0]
    output:
        demuxr1=join(config["demuxdir"],"{sample}-R1.fastq.gz"),
        demuxr2=join(config["demuxdir"],"{sample}-R2.fastq.gz"),
        demuxlog=join(config["demuxdir"],"{sample}.extract.log")
    conda:
        "../../workflow/envs/demultiplexBarcodesNoBuilds.yml"
    shell:
        """
        umi_tools extract --extract-method string \
        --bc-pattern NNNNNNN --bc-pattern2 NNNNNNN -L {output.demuxlog} \
        -I {input.r1} --read2-in {input.r2} \
        -S {output.demuxr1} --read2-out {output.demuxr2}
        """

# Use grep to split out reads that match each known barcode combination.
rule splitReads:
    input:
        demuxr1=join(config["demuxdir"],"{sample}-R1.fastq.gz"),
        demuxr2=join(config["demuxdir"],"{sample}-R2.fastq.gz"),
        indices=config["indices"]
    params:
        demuxdir=config["demuxdir"],
        demux=config["demux"]
    output:
        expand(join(config["demuxdir"],"R1/{{sample}}-L{phase}.fastq.gz"),
            phase=phases, allow_missing=True),
        expand(join(config["demuxdir"],"R2/{{sample}}-L{phase}.fastq.gz"),
            phase=phases, allow_missing=True)
    shell:
        """
        {params.demux} -x {input.demuxr1} -y {input.demuxr2} \
          -o {params.demuxdir} -s {wildcards.sample} -i {input.indices}
        """

# Use cutadapt to trim the remaining primer sequences from the sequencing reads.
# Note that the length of these sequences depends on the round 1 index that was used.
# Also rename the files to match the sample name.
rule trimReads:
    input:
        demuxr1="".join([config["demuxdir"],"/R1/{file}.fastq.gz"]),
        demuxr2="".join([config["demuxdir"],"/R2/{file}.fastq.gz"])
    params:
        trimdir = join(config["trimdir"],"{group}"),
        trim1= lambda wildcards: str(int(config["lenR1primer"])-int(config["lenR1index"])+int((wildcards.file).split("-L")[1])),
        trim2= lambda wildcards: str(int(config["lenR2primer"])-int((wildcards.file).split("-L")[1]))
    output:
        trimr1= join(config["trimdir"], "{group}/R1/{file}_{sample}.fastq.gz"),
        trimr2= join(config["trimdir"], "{group}/R2/{file}_{sample}.fastq.gz")
    conda:
        "../../workflow/envs/demultiplexBarcodesNoBuilds.yml"
    shell:
        """
        mkdir -p {params.trimdir}
        mkdir -p {params.trimdir}/R1
        mkdir -p {params.trimdir}/R2
        cutadapt -u {params.trim1} -U {params.trim2} \
          -o {output.trimr1} -p {output.trimr2} \
          {input.demuxr1} {input.demuxr2}
        """

# Count the number of sequencing reads in each file.
rule countReads:
	params:
		dir="".join([config["trimdir"],"/{group}/R1"])
	output:
		readSummary="".join([config["trimdir"],"/{group}/summary.txt"])
	shell:
		"""
		rm -f {output.readSummary}
		for f in {params.dir}/*.fastq.gz
		do
			numReads="$( zcat ${{f}} | grep ^@ | wc -l )"
			echo ${{f}} ${{numReads}}
		done > {output.readSummary}
		"""

# Remove files with low read numbers.
rule removeLowReads:
    params:
        readSummary="".join([config["trimdir"],"/{group}/summary.txt"]),
        lowReadsSummary="".join([config["trimdir"],"/{group}/lowReadsSummary.txt"]),
        removedDir="".join([config["trimdir"],"/{group}/removed"])
    output:
        removedFiles="".join([config["trimdir"],"/{group}/lowReadsSummary.txt"])
    shell:
        """
        less {params.readSummary} | sort -n -k2 | awk '$2 <100' | cut -f1 -d' ' > {params.lowReadsSummary}
        mkdir -p {params.removedDir}
        while read line
        do
            if [ -f $line ]; then
                mv $line {params.removedDir}
                line=${{line/R1/R2}}
                mv $line {params.removedDir}
            fi
        done < {params.lowReadsSummary}
        """
