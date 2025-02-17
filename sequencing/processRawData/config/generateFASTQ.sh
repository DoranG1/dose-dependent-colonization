paste <( ls /oak/stanford/groups/relman/raw_data/dorang/20210819-e0015-e0020-16S/16S_*R1* ) <( ls /oak/stanford/groups/relman/raw_data/dorang/20210819-e0015-e0020-16S/16S_*R2* ) \
    > fastq-raw.txt
while read r1 r2
do
    sample=${r1##*/}
    sample=${sample%_*-*.fastq.gz}
    sample=${sample/_/-}
    sample=${sample/_/-}
    echo ${r1} ${r2} ${sample}
done < fastq-raw.txt > fastq.txt
rm -f fastq-raw.txt
